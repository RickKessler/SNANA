// Created Apr 2025
// Functions to read/write npz file used by python.
//
// Oct 28 2026:
//    Fix string columns (CID, FIELD): header was written as int8 2D array
//    (cnpy maps char -> '|i1'), so python read back integers. Now write
//    1D fixed-width byte-string header, dtype '|S<wsize>'.
//
// 
#define MXTABLE_NPZ  10
#define ISFILE_NPZ_COV        1
#define ISFILE_NPZ_FACTORIZE  2


#ifdef __cplusplus
extern"C" {
#endif

  void INIT_NPZFILES(char *PREFIX);
  void OPEN_NPZFILE(char *FILENAME);
  void SNTABLE_CREATE_NPZ(int IDTABLE, char *TBNAME, char *PGNAME);
  void SNTABLE_ADDCOL_NPZ(int IDTABLE, void *PTRVAR, 
			  SNTABLE_ADDCOL_VARDEF *ADDCOL_VARDEF) ;

  void SNTABLE_FILL_NPZ(int IDTABLE);
  void CLOSE_NPZFILE(void);

  // Nov 2026: read table from npz file (analog of SNTABLE_xxx_TEXT)
  int  SNTABLE_NEVT_NPZ(char *FILENAME);
  int  SNTABLE_READPREP_NPZ(void);
  int  SNTABLE_READ_EXEC_NPZ(int OPTMASK);

  int  get_ITABLE_SPARSE(int IDTABLE, char *FUNNAM, int OPT_ABORT );

  int read_npz_covmat(char *npz_file, double *array1d);
  int read_npz_factorized(char *npz_file, double **diag, double **U, int *K);
  int get_npz_filetype(char *npz_file);

  // ---- table writer (buffered, spilled to compressed temp files, zipped at close)
  int  npz_table_create (char *fileName);
  int  npz_table_addcol (int h, char *varName, int icast, int strlen_max);
  void npz_table_append (int h, int icol, void *ptr);
  void npz_table_endrow (int h);     // call once per row; triggers flush every NFLUSH rows
  void npz_table_close  (int h);
  void npz_table_set_header(int h, const char *json); // JSON string -> zip archive comment

  void  errmsg ( int isev, int iprompt, char *fnam, char *msg1, char *msg2 );
  void  print_banner ( const char *banner ) ;
  void  print_preAbort_banner(char *fnam);
  void  trim_blank_spaces(char *string) ;
  int   strcmp_ignoreCase(char *str1, char *str2) ;
  void  debugexit(char *string);
  void catVarList_with_comma(char *varList, char *addVarName);
  int  IGNOREFILE(char *fileName);

#ifdef __cplusplus
}          
#endif
// Created Apr 2025
// Utilities to read/write python-style npz files 
// Initial use is for wfit.c to read npz-formatted covariance file 
// from create_covariance.py.
// Use package from https://github.com/rogersce/cnpy


// ============================================
//           DECLARATIONS
// ============================================

int  NTABLE_NPZ ;
char FILEPREFIX_NPZ[MXPATHLEN];

#define SUFFIX_NPZ   "NPZ"   
#define MXTABLE_NPZ 10

struct TABLEINFO_NPZ {
  int    NTABLE ;
  int    IDTABLE[MXTABLE_NPZ] ;
  char   TBNAME[MXTABLE_NPZ][40] ;  // table name; e.g., SNANA or FITRES or LCPLOT
  char   PGNAME[MXTABLE_NPZ][300] ; 

  char   FILENAME[MXTABLE_NPZ][MXCHAR_FILENAME] ;
  int    NFILL[MXTABLE_NPZ] ; // number of fill calls per table

  int    HANDLE[MXTABLE_NPZ] ;  // handle returned by npz_table_create

  // per-column info, filled by SNTABLE_ADDCOL_NPZ (analog of TABLEINFO_TEXT)
  int     NVAR[MXTABLE_NPZ];
  char  **VARNAME[MXTABLE_NPZ];  // [ITAB][IVAR] column name
  int    *ICAST[MXTABLE_NPZ] ;   // [ITAB][IVAR] ICAST_[D,F,I,S,L,C]
  int    *ICOL[MXTABLE_NPZ] ;    // [ITAB][IVAR] column index from npz_table_addcol
  void  **PTRVAR[MXTABLE_NPZ] ;  // [ITAB][IVAR] address of variable to write

} TABLEINFO_NPZ ;

// ====================================
void INIT_NPZFILES(char *PREFIX) {

  // Oct 2026
  // store prefix, and misc inits.
  // Cannot open table file(s) here because we don't know
  // which one(s). Table File(s) opened at SNTABLE_CREATE stage.

  int LEN_PREFIX = strlen(PREFIX);
  //  char *FMT, comment[200];

  char fnam[] = "INIT_NPZFILES" ;  (void)fnam;

  // --------- BEGIN -----------

  if ( LEN_PREFIX >= MXPATHLEN ) {
    print_preAbort_banner(fnam);
    printf("  PREFIX = '%s' \n", PREFIX);
    sprintf(MSGERR1, "LEN(PREFIX) = %d exceeds bound MXPATHLEN=%d", LEN_PREFIX, MXPATHLEN);
    sprintf(MSGERR2, "Cannot copy to FILEPREFIX_NPZ");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  sprintf(FILEPREFIX_NPZ, "%s", PREFIX);

  TABLEINFO_NPZ.NTABLE = 0 ;
 
  return;

} // end INIT_NPZFILES

// =====================================================================
//   Nov 2026: NPZ table READER (analog of the _TEXT read functions).
//
//   Reads the npz table written by SNTABLE_xxx_NPZ (or by python
//   numpy.savez/savez_compressed): one 1-D .npy per column, column order
//   taken from the zip central directory, header dictionary taken from
//   the zip comment (JSON string:string). Columns are inflated on demand
//   by SNTABLE_READ_EXEC_NPZ. Zip64 (>4GB) tables are not supported.
// =====================================================================
#include <vector>
#include <string>
#include <algorithm>
#include <strings.h>
#include <stdint.h>
#include <zlib.h>

namespace npzrd {

  struct COL {
    std::string name;
    int      icast;       // ICAST_[D,F,I,S,L,C] of data stored in file
    int      kind;        // 'f' 'i' 'u' 'S' 'U'
    int      wsize;       // bytes per element in file
    size_t   nrow;
    uint32_t method, csize, usize;
    uint64_t data_off;    // file offset of (compressed) data stream
    size_t   npy_off;     // offset of first element inside inflated stream
  };

  struct INFO {
    std::string file;
    std::vector<COL> col;
    std::string header_json;
    bool opened;
    INFO() : opened(false) {}
  } R;

  uint16_t u16(const unsigned char *p) { return (uint16_t)(p[0] | (p[1]<<8)); }
  uint32_t u32(const unsigned char *p) 
  { return (uint32_t)p[0] | ((uint32_t)p[1]<<8) | ((uint32_t)p[2]<<16) | ((uint32_t)p[3]<<24); }

  // inflate (at most maxout bytes of) the column stream; maxout=0 -> all
  bool inflate_col(FILE *fp, const COL &c, std::vector<char> &out, size_t maxout) {
    size_t want = c.usize ;
    if ( maxout > 0 && maxout < want ) { want = maxout; }
    out.assign(want, 0);
    if ( want == 0 ) { return true; }
    if ( fseeko(fp, (off_t)c.data_off, SEEK_SET) != 0 ) { return false; }

    if ( c.method == 0 ) { return fread(&out[0], 1, want, fp) == want; }
    if ( c.method != 8 ) { return false; }

    z_stream zs;  memset(&zs, 0, sizeof(zs));
    if ( inflateInit2(&zs, -15) != Z_OK ) { return false; }
    std::vector<unsigned char> in(65536);
    size_t left = c.csize;
    zs.next_out  = (Bytef*)&out[0];
    zs.avail_out = (uInt)want;
    int ret = Z_OK;
    while ( zs.avail_out > 0 && ret == Z_OK ) {
      if ( zs.avail_in == 0 ) {
	size_t n = std::min(left, in.size());
	if ( n == 0 ) { break; }
	if ( fread(&in[0], 1, n, fp) != n ) { inflateEnd(&zs); return false; }
	left -= n;  zs.next_in = &in[0];  zs.avail_in = (uInt)n;
      }
      ret = inflate(&zs, Z_NO_FLUSH);
    }
    bool ok = ( zs.avail_out == 0 || ret == Z_STREAM_END );
    inflateEnd(&zs);
    return ok;
  }

  // parse the npy header at start of inflated stream
  bool parse_npy(const std::vector<char> &b, COL &c, std::string &err) {
    const unsigned char *p = (const unsigned char*)&b[0];
    if ( b.size() < 12 || p[0] != 0x93 || memcmp(p+1,"NUMPY",5) != 0 )
      { err = "bad npy magic in column " + c.name; return false; }
    int major = p[6];
    size_t hlen, start;
    if ( major == 1 ) { hlen = u16(p+8); start = 10; }
    else              { hlen = u32(p+8); start = 12; }
    if ( start + hlen > b.size() ) { err = "npy header too long: " + c.name; return false; }
    std::string h(b.begin()+start, b.begin()+start+hlen);
    c.npy_off = start + hlen;

    size_t k = h.find("descr");
    if ( k == std::string::npos ) { err = "no descr in " + c.name; return false; }
    k = h.find(':', k);  k = h.find('\'', k);
    if ( k == std::string::npos ) { err = "bad descr in " + c.name; return false; }
    size_t k2 = h.find('\'', k+1);
    std::string d = h.substr(k+1, k2-k-1);        // e.g. <f8  |S20  <i4
    if ( d.size() < 3 || d[0] == '>' ) 
      { err = "unsupported dtype '" + d + "' in column " + c.name; return false; }
    c.kind  = d[1];
    c.wsize = atoi(d.c_str()+2);

    k = h.find("fortran_order");
    k = h.find('(', k == std::string::npos ? 0 : k);
    if ( k == std::string::npos ) { err = "no shape in " + c.name; return false; }
    size_t k3 = h.find(')', k);
    std::string shp = h.substr(k+1, k3-k-1);
    if ( shp.find(',') != std::string::npos && shp.find_first_of("0123456789",shp.find(',')) != std::string::npos )
      { err = "column " + c.name + " is not 1-D"; return false; }
    c.nrow = (size_t) atoll(shp.c_str());

    switch ( c.kind ) {
    case 'f': c.icast = (c.wsize==8) ? ICAST_D : ICAST_F ;  break;
    case 'i': c.icast = (c.wsize==8) ? ICAST_L : (c.wsize==4) ? ICAST_I : ICAST_S ; break;
    case 'u': c.icast = (c.wsize>=4) ? ICAST_L : (c.wsize==2) ? ICAST_I : ICAST_S ; break;
    case 'S': case 'U': c.icast = ICAST_C ; break;
    default:  err = std::string("unsupported dtype '") + d + "' in column " + c.name; return false;
    }
    if ( c.kind=='f' && c.wsize!=8 && c.wsize!=4 ) 
      { err = "unsupported float size in " + c.name; return false; }
    return true;
  }

  bool open_file(const char *fileName, std::string &err) {
    R = INFO();  R.file = fileName;
    FILE *fp = fopen(fileName, "rb");
    if ( !fp ) { err = "cannot open file"; return false; }

    // ---- end-of-central-directory record (search back for signature)
    fseeko(fp, 0, SEEK_END);
    off_t fsize = ftello(fp);
    size_t ntail = (size_t) std::min<off_t>(fsize, 22 + 65535);
    std::vector<unsigned char> tail(ntail);
    fseeko(fp, fsize - (off_t)ntail, SEEK_SET);
    if ( fread(&tail[0], 1, ntail, fp) != ntail ) { fclose(fp); err = "read error"; return false; }
    int ie = -1;
    for ( int i = (int)ntail - 22; i >= 0; i-- ) {
      if ( tail[i]==0x50 && tail[i+1]==0x4b && tail[i+2]==0x05 && tail[i+3]==0x06 ) { ie = i; break; }
    }
    if ( ie < 0 ) { fclose(fp); err = "not a zip/npz file"; return false; }
    const unsigned char *e = &tail[ie];
    uint16_t ncd = u16(e+10);   uint32_t cdoff = u32(e+16);  uint16_t lcom = u16(e+20);
    if ( ncd == 0xFFFF || cdoff == 0xFFFFFFFFu ) { fclose(fp); err = "zip64 npz not supported"; return false; }
    if ( lcom > 0 && ie + 22 + lcom <= (int)ntail )
      { R.header_json.assign((const char*)e+22, lcom); }

    // ---- central directory (gives column order)
    fseeko(fp, (off_t)cdoff, SEEK_SET);
    std::vector<uint32_t> lho;
    for ( int i=0; i < ncd; i++ ) {
      unsigned char ch[46];
      if ( fread(ch,1,46,fp) != 46 || u32(ch) != 0x02014b50u ) { fclose(fp); err = "bad central directory"; return false; }
      COL c;
      c.method = u16(ch+10);  c.csize = u32(ch+20);  c.usize = u32(ch+24);
      uint16_t nl = u16(ch+28), xl = u16(ch+30), cl = u16(ch+32);
      uint32_t lo = u32(ch+42);
      if ( c.csize==0xFFFFFFFFu || c.usize==0xFFFFFFFFu || lo==0xFFFFFFFFu )
	{ fclose(fp); err = "zip64 npz not supported"; return false; }
      std::string nm(nl,' ');
      if ( fread(&nm[0],1,nl,fp) != nl ) { fclose(fp); err = "read error"; return false; }
      fseeko(fp, xl + cl, SEEK_CUR);
      if ( nm.size() > 4 && nm.compare(nm.size()-4,4,".npy") == 0 ) { nm.erase(nm.size()-4); }
      c.name = nm;  c.icast = -9;  c.kind = 0;  c.wsize = 0;  c.nrow = 0;  c.data_off = 0;  c.npy_off = 0;
      R.col.push_back(c);  lho.push_back(lo);
    }

    // ---- local headers -> data offset; npy header -> dtype, nrow
    for ( size_t i=0; i < R.col.size(); i++ ) {
      unsigned char lh[30];
      fseeko(fp, (off_t)lho[i], SEEK_SET);
      if ( fread(lh,1,30,fp) != 30 || u32(lh) != 0x04034b50u ) { fclose(fp); err = "bad local header"; return false; }
      R.col[i].data_off = (uint64_t)lho[i] + 30 + u16(lh+26) + u16(lh+28);
      std::vector<char> b;
      if ( !inflate_col(fp, R.col[i], b, 1024) ) { fclose(fp); err = "cannot inflate column " + R.col[i].name; return false; }
      if ( !parse_npy(b, R.col[i], err) ) { fclose(fp); return false; }
    }
    fclose(fp);

    if ( R.col.empty() ) { err = "no columns in npz file"; return false; }
    for ( size_t i=1; i < R.col.size(); i++ ) 
      if ( R.col[i].nrow != R.col[0].nrow ) { err = "columns have different lengths: " + R.col[i].name; return false; }

    // make sure the ID column is first (CID, CCID, GALID, ...)
    {
      const char *idlist[] = VARNAME_TABLE_IDLIST ;
      int found = -1;
      for ( int j=0; j < N_VARNAME_TABLE_IDLIST && found < 0; j++ )
	for ( size_t i=0; i < R.col.size(); i++ )
	  if ( strcasecmp(R.col[i].name.c_str(), idlist[j]) == 0 ) { found = (int)i; break; }   // case-insensitive: SNID, snid, CID ...
      if ( found > 0 ) { COL t = R.col[found]; R.col.erase(R.col.begin()+found); R.col.insert(R.col.begin(), t); }
    }
    R.opened = true;
    return true;
  }

  // parse flat JSON {"key": "value", ...} from zip comment -> "# KEY: value"
  void json_to_comments(std::vector<std::string> &lines) {
    const std::string &j = R.header_json;
    size_t i = 0, n = j.size();
    std::vector<std::string> items;
    while ( i < n ) {
      if ( j[i] != '"' ) { i++; continue; }
      std::string s;  i++;
      while ( i < n && j[i] != '"' ) {
	if ( j[i] == '\\' && i+1 < n ) {
	  i++;
	  switch ( j[i] ) {
	  case 'n': s += ' '; break;  case 't': s += ' '; break;  case 'r': break;
	  case 'u': if ( i+4 < n ) { s += (char) strtol(j.substr(i+1,4).c_str(), NULL, 16); i += 4; } break;
	  default:  s += j[i];
	  }
	} else { s += j[i]; }
	i++;
      }
      i++;  items.push_back(s);
    }
    for ( size_t k=0; k+1 < items.size(); k += 2 ) {
      std::string key = items[k];
      if ( key.find(':') == std::string::npos ) { key += ":"; }
      lines.push_back("# " + key + " " + items[k+1] + "\n");
    }
  }

} // end namespace npzrd

static bool npzrd_open(char *fileName, std::string &err) { return npzrd::open_file(fileName, err); }

// ==============================================
int SNTABLE_NEVT_NPZ(char *FILENAME) {
  // return number of rows; open file if FILENAME is non-blank, else use opened file
  std::string err;
  if ( strlen(FILENAME) > 0 && !( npzrd::R.opened && npzrd::R.file == FILENAME ) ) {
    if ( !npzrd::open_file(FILENAME, err) ) { return 0; }
  }
  if ( !npzrd::R.opened ) { return 0; }
  return (int) npzrd::R.col[0].nrow ;
} // end SNTABLE_NEVT_NPZ

// ==============================================
int SNTABLE_READPREP_NPZ(void) {

  // Analog of SNTABLE_READPREP_TEXT: fill VARNAME and ICAST for each column.
  // ID column (first) is always stored as string, even if integers in file.

  int NVAR = (int) npzrd::R.col.size();
  char fnam[] = "SNTABLE_READPREP_NPZ" ;

  if ( !npzrd::R.opened ) {
    sprintf(MSGERR1, "npz file is not open");  sprintf(MSGERR2, "Call TABLEFILE_OPEN first");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }
  if ( NVAR >= MXVAR_TABLE ) {
    sprintf(MSGERR1, "NVAR=%d exceeds bound MXVAR_TABLE=%d", NVAR, MXVAR_TABLE);
    sprintf(MSGERR2, "%.200s", npzrd::R.file.c_str());
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  for ( int ivar=0; ivar < NVAR; ivar++ ) {
    npzrd::COL &c = npzrd::R.col[ivar];
    snprintf(READTABLE_POINTERS.VARNAME[ivar], MXCHAR_VARNAME, "%s", c.name.c_str());
    int ICAST = c.icast;
    if ( ivar == 0 ) { ICAST = ICAST_C; }
    READTABLE_POINTERS.ICAST_READ[ivar]  = ICAST ;
    READTABLE_POINTERS.ICAST_STORE[ivar] = ICAST ;
  }
  READTABLE_POINTERS.NVAR_TOT = NVAR ;
  return NVAR ;

} // end SNTABLE_READPREP_NPZ

// ==============================================
int SNTABLE_READ_EXEC_NPZ(int OPTMASK) {

  // Analog of SNTABLE_READ_EXEC_TEXT: inflate each requested column and fill
  // the user arrays set by SNTABLE_READPREP_VARDEF. Returns NROW.
  // OPTMASK bit0: store header comments (from zip comment).

  bool STORE_TABLECOMMENT = (OPTMASK & 1) > 0 ;
  int  NVAR_READ = READTABLE_POINTERS.NVAR_READ ;
  int  NROW      = (int) npzrd::R.col[0].nrow ;
  char fnam[]    = "SNTABLE_READ_EXEC_NPZ" ;

  if ( NROW > READTABLE_POINTERS.MXLEN && READTABLE_POINTERS.MXLEN > 0 ) {
    sprintf(MSGERR1, "NROW=%d exceeds memory bound MXLEN=%d", NROW, READTABLE_POINTERS.MXLEN);
    sprintf(MSGERR2, "%.200s", npzrd::R.file.c_str());
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  if ( STORE_TABLECOMMENT ) {
    std::vector<std::string> lines;
    npzrd::json_to_comments(lines);
    READTABLE_POINTERS.NLINE_TABLECOMMENT = 0 ;
    for ( size_t k=0; k < lines.size() && k < (size_t)MXLINE_TABLECOMMENT; k++ ) 
      { store_read_comment_TEXT((char*)lines[k].c_str()); }  // same storage as text
  }

  FILE *fp = fopen(npzrd::R.file.c_str(), "rb");
  if ( !fp ) {
    sprintf(MSGERR1, "Cannot open %.200s", npzrd::R.file.c_str());  MSGERR2[0] = 0;
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  std::vector<char> buf;
  for ( int i=0; i < NVAR_READ; i++ ) {
    int ivar = READTABLE_POINTERS.PTRINDEX[i];
    npzrd::COL &c = npzrd::R.col[ivar];
    int ICAST = READTABLE_POINTERS.ICAST_STORE[ivar];
    int NPTR  = READTABLE_POINTERS.NPTR[ivar];
    if ( NPTR <= 0 || READTABLE_POINTERS.MXLEN == 0 ) { continue; }

    if ( !npzrd::inflate_col(fp, c, buf, 0) || buf.size() < c.npy_off + (size_t)c.wsize*c.nrow ) {
      sprintf(MSGERR1, "Could not read column '%s'", c.name.c_str());
      sprintf(MSGERR2, "%.200s", npzrd::R.file.c_str());
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
    }
    const char *base = &buf[0] + c.npy_off;

    for ( int isn=0; isn < NROW; isn++ ) {
      const char *p = base + (size_t)isn * c.wsize;
      long double dv = 0.0L ;  char cv[MXCHAR_VARNAME+8];  cv[0] = 0;
      
      if ( c.kind == 'S' || c.kind == 'U' ) {
	// string: take low byte of each UCS4 char for 'U'
	int nch = (c.kind=='S') ? c.wsize : c.wsize/4, n=0;
	for ( int j=0; j < nch && n < MXCHAR_VARNAME-1; j++ ) {
	  char ch = (c.kind=='S') ? p[j] : p[4*j];
	  if ( ch == 0 ) { break; }
	  cv[n++] = ch;
	}
	while ( n > 0 && (cv[n-1]==' ' || cv[n-1]=='\t') ) { n--; }  // strip trailing blanks (as text parsing does)
	cv[n] = 0;  dv = strtold(cv, NULL);
      }
      else {
	if      ( c.kind=='f' && c.wsize==8 ) { double x; memcpy(&x,p,8); dv = x; }
	else if ( c.kind=='f' )               { float  x; memcpy(&x,p,4); dv = x; }
	else if ( c.wsize==8 )  { long long x; memcpy(&x,p,8); dv = (long double)x; }
	else if ( c.wsize==4 )  { if ( c.kind=='i' ) {int x; memcpy(&x,p,4); dv = x;} else {unsigned x; memcpy(&x,p,4); dv = x;} }
	else if ( c.wsize==2 )  { if ( c.kind=='i' ) {short x; memcpy(&x,p,2); dv = x;} else {unsigned short x; memcpy(&x,p,2); dv = x;} }
	else                    { dv = (c.kind=='i') ? (long double)(signed char)p[0] : (long double)(unsigned char)p[0]; }
	if ( ICAST == ICAST_C ) { sprintf(cv, "%lld", (long long)dv); }   // numeric ID -> string
      }

      for ( int nptr=0; nptr < NPTR; nptr++ ) {
	if      ( ICAST == ICAST_D ) { READTABLE_POINTERS.PTRVAL_D[nptr][ivar][isn] = (double)dv; }
	else if ( ICAST == ICAST_F ) { READTABLE_POINTERS.PTRVAL_F[nptr][ivar][isn] = (float)dv; }
	else if ( ICAST == ICAST_I ) { READTABLE_POINTERS.PTRVAL_I[nptr][ivar][isn] = (int)dv; }
	else if ( ICAST == ICAST_S ) { READTABLE_POINTERS.PTRVAL_S[nptr][ivar][isn] = (short int)dv; }
	else if ( ICAST == ICAST_L ) { READTABLE_POINTERS.PTRVAL_L[nptr][ivar][isn] = (long long int)dv; }
	else if ( ICAST == ICAST_C ) { sprintf(READTABLE_POINTERS.PTRVAL_C[nptr][ivar][isn], "%s", cv); }
	else {
	  sprintf(MSGERR1,"Unknown ICAST=%d  var=%s", ICAST, READTABLE_POINTERS.VARNAME[ivar]);
	  sprintf(MSGERR2,"See ICAST_  parameters in sntools_output.h");
	  errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2 );
	}
      }
    } // isn
  } // i
  fclose(fp);

  // allow opening another file for reading
  NAME_TABLEFILE[OPENFLAG_READ][IFILETYPE_NPZ][0] = 0 ;
  USE_TABLEFILE[OPENFLAG_READ][IFILETYPE_NPZ]     = 0 ;
  npzrd::R.opened = false;

  return NROW ;

} // end SNTABLE_READ_EXEC_NPZ

// ==============================================
void OPEN_NPZFILE(char *FILENAME) {

  // Nov 2026: open npz table file for reading. Parse zip directory and
  // npy headers (column names, types, nrow) and zip comment (header dict).
  // Column data are inflated later by SNTABLE_READ_EXEC_NPZ.

  char fnam[] = "OPEN_NPZFILE" ; 
  std::string err;

  // ----------- BEGIN --------

  if ( !npzrd_open(FILENAME, err) ) {
    print_preAbort_banner(fnam);
    printf("  file = %s\n", FILENAME);
    sprintf(MSGERR1, "Could not read npz table file:");
    sprintf(MSGERR2, "%.200s", err.c_str());
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  return ;
} // end OPEN_NPZFILE


// ==============================================
void SNTABLE_CREATE_NPZ(int IDTABLE, char *TBNAME, char *PGNAME) {

  // Created Oct 2026
  // IDTABLE = int table id
  // NAME    = char table name (e.g., SNANA or FITRES) used to construct filename
  // PGNAME  = name of calling program to be written to data header
  //
  // Analog of SNTABLE_CREATE_TEXT. File name is [PREFIX].[TBNAME].NPZ.
  // Register the table and allocate per-column storage; columns are
  // defined later by SNTABLE_ADDCOL_NPZ. Unlike the text table, the file
  // is not written until the table is closed (npz_table_close), so there
  // is no file pointer to store here.

  int  ITAB = TABLEINFO_NPZ.NTABLE ;
  int  ivar, MX = MXVAR_TABLE ;
  char FILENAME[MXPATHLEN];
  char fnam[] = "SNTABLE_CREATE_NPZ" ;  (void)fnam;

  // ------------- BEGIN --------------

  // table arrays are fixed size; abort before indexing past them
  if ( ITAB >= MXTABLE_NPZ ) {
    sprintf(MSGERR1, "Cannot create table %s: too many tables", TBNAME);
    sprintf(MSGERR2, "MXTABLE_NPZ = %d", MXTABLE_NPZ);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  printf("  %s: init %s %s-table (ITAB=%d). \n", fnam, TBNAME, SUFFIX_NPZ, ITAB); 
  fflush(stdout);

  // construct file name from stored prefix and table name.
  if ( strchr(FILEPREFIX_NPZ,'.') == NULL ) {   
    // no dot --> add on table and suffix to file name
    sprintf(FILENAME, "%.*s.%s.%s", MXPATHLEN-20, FILEPREFIX_NPZ, TBNAME, SUFFIX_NPZ);
  }
  else {
    // if there is a dot in the prefix, assume it's the full filename
    sprintf(FILENAME, "%s", FILEPREFIX_NPZ);
  }

  // store other info (same as for TEXT format)
  TABLEINFO_NPZ.IDTABLE[ITAB] = IDTABLE ;
  TABLEINFO_NPZ.NFILL[ITAB]   = 0 ;
  TABLEINFO_NPZ.NVAR[ITAB]    = 0 ;
  sprintf(TABLEINFO_NPZ.TBNAME[ITAB],   "%s", TBNAME);
  sprintf(TABLEINFO_NPZ.PGNAME[ITAB],   "%s", PGNAME);
  sprintf(TABLEINFO_NPZ.FILENAME[ITAB], "%s", FILENAME);
  TABLEINFO_NPZ.NTABLE++ ;

  // open table (buffered writer; returns table handle)
  int itab = npz_table_create(FILENAME); 
  // abort if returned itab index does not match expected ITAB index
  if ( itab != ITAB ) {
    sprintf(MSGERR1, "Expect ITAB = %d", ITAB);
    sprintf(MSGERR2, "but npz_table_create returns itab=%d", itab);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }
  TABLEINFO_NPZ.HANDLE[ITAB] = itab ;

  // allocate memory for column names, casts, column indices and pointers
  TABLEINFO_NPZ.VARNAME[ITAB] = (char**) malloc( MX * sizeof(char*) );
  for(ivar=0; ivar < MX; ivar++ ) 
    { TABLEINFO_NPZ.VARNAME[ITAB][ivar] = (char*)malloc(MXCHAR_VARNAME*sizeof(char)); }

  TABLEINFO_NPZ.ICAST[ITAB]  = (int*)   malloc( MX * sizeof(int)   );
  TABLEINFO_NPZ.ICOL[ITAB]   = (int*)   malloc( MX * sizeof(int)   );
  TABLEINFO_NPZ.PTRVAR[ITAB] = (void**) malloc( MX * sizeof(void*) );

  for(ivar=0; ivar < MX; ivar++ ) { 
    TABLEINFO_NPZ.ICAST[ITAB][ivar]  = -9 ; 
    TABLEINFO_NPZ.ICOL[ITAB][ivar]   = -9 ; 
    TABLEINFO_NPZ.PTRVAR[ITAB][ivar] = NULL ; 
  }

  return ;

} // end SNTABLE_CREATE_NPZ


// ==============================================
void SNTABLE_ADDCOL_NPZ(int IDTABLE, void *PTRVAR, 
			SNTABLE_ADDCOL_VARDEF *ADDCOL_VARDEF) {

  // Analog of SNTABLE_ADDCOL_TEXT.
  // For each variable in ADDCOL_VARDEF, define a column in the npz table
  // (npz_table_addcol) and store its name, cast, column index, and address.
  // The address is read at each SNTABLE_FILL_NPZ call.
  // All columns must be added before the first FILL.

  int  ITAB, IVAR, ivar, ICAST, ICOL, HANDLE, STRLEN ;
  char VARNAME[MXCHAR_VARNAME], *TBNAME ;
  char fnam[] = "SNTABLE_ADDCOL_NPZ";

  // ----------- BEGIN -----------

  // find sparse index ITAB for this IDTABLE
  ITAB   = get_ITABLE_SPARSE(IDTABLE, fnam, 0 );

  if ( ITAB < 0 ) {
    print_preAbort_banner(fnam);
    printf("   VARLIST_ORIG = '%s' \n", ADDCOL_VARDEF->VARLIST_ORIG);
    printf("   NVAR = %d \n", ADDCOL_VARDEF->NVAR );
    sprintf(MSGERR1, "Could not find npz table with ");
    sprintf(MSGERR2, "IDTABLE = %d", IDTABLE ); 
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
  }

  TBNAME = TABLEINFO_NPZ.TBNAME[ITAB] ;
  HANDLE = TABLEINFO_NPZ.HANDLE[ITAB] ;

  if ( TABLEINFO_NPZ.NFILL[ITAB] > 0 ) {
    sprintf(MSGERR1, "Cannot add columns to %s table after rows are filled", TBNAME);
    sprintf(MSGERR2, "VARLIST_ORIG = '%.*s'", MXCHAR_MSGERR-20, ADDCOL_VARDEF->VARLIST_ORIG);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  for(ivar=0 ; ivar < ADDCOL_VARDEF->NVAR; ivar++ ) {

    IVAR = TABLEINFO_NPZ.NVAR[ITAB] ;   // global column counter for this table

    // unlike TEXT (which silently skips), abort since npz columns must be consistent
    if ( IVAR >= MXVAR_TABLE ) {
      sprintf(MSGERR1, "NVAR exceeds bound of MXVAR_TABLE=%d for %s table", 
	      MXVAR_TABLE, TBNAME);
      sprintf(MSGERR2, "while adding VARNAME='%s'", ADDCOL_VARDEF->VARNAME[ivar]);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
    }

    sprintf(VARNAME, "%s", ADDCOL_VARDEF->VARNAME[ivar] );    

    // replace CCID -> CID, same as for text file
    if ( strcmp(VARNAME,"CCID") == 0 )  { sprintf(VARNAME,"CID") ; }

    ICAST = ADDCOL_VARDEF->ICAST[ivar] ;

    // string width in bytes (only used for ICAST_C)
    STRLEN = 0 ;
    if ( ICAST == ICAST_C ) {
      STRLEN = ADDCOL_VARDEF->ISIZE[ivar] ;
      if ( STRLEN <= 0 ) { STRLEN = MXCHAR_CCID ; }
    }

    // define column in the npz writer; returns column index
    ICOL = npz_table_addcol(HANDLE, VARNAME, ICAST, STRLEN);
    if ( ICOL < 0 ) {
      sprintf(MSGERR1, "Invalid ICAST=%d for VARNAME='%s' (%s table)", 
	      ICAST, VARNAME, TBNAME);
      sprintf(MSGERR2, "Valid ICAST: D=%d F=%d I=%d S=%d L=%d C=%d", 
	      ICAST_D, ICAST_F, ICAST_I, ICAST_S, ICAST_L, ICAST_C);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
    }

    // npz column index must track our own counter; FILL relies on this
    if ( ICOL != IVAR ) {
      sprintf(MSGERR1, "Expect ICOL=%d for VARNAME='%s'", IVAR, VARNAME);
      sprintf(MSGERR2, "but npz_table_addcol returns ICOL=%d", ICOL);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
    }

    // store column info
    sprintf(TABLEINFO_NPZ.VARNAME[ITAB][IVAR], "%s", VARNAME);
    TABLEINFO_NPZ.ICAST[ITAB][IVAR] = ICAST ;
    TABLEINFO_NPZ.ICOL[ITAB][IVAR]  = ICOL ;

    // store pointer based on cast (same convention as SNTABLE_ADDCOL_TEXT):
    // numeric PTRVAR is an array indexed by ivar; string PTRVAR is one char array.
    if ( ICAST == ICAST_D ) 
      { TABLEINFO_NPZ.PTRVAR[ITAB][IVAR] = (void*)((double*)PTRVAR + ivar) ; }
    else if ( ICAST == ICAST_F ) 
      { TABLEINFO_NPZ.PTRVAR[ITAB][IVAR] = (void*)((float*)PTRVAR + ivar) ; }
    else if ( ICAST == ICAST_I ) 
      { TABLEINFO_NPZ.PTRVAR[ITAB][IVAR] = (void*)((int*)PTRVAR + ivar) ; }
    else if ( ICAST == ICAST_S ) 
      { TABLEINFO_NPZ.PTRVAR[ITAB][IVAR] = (void*)((short int*)PTRVAR + ivar) ; }
    else if ( ICAST == ICAST_L ) 
      { TABLEINFO_NPZ.PTRVAR[ITAB][IVAR] = (void*)((long long int*)PTRVAR + ivar) ; }
    else if ( ICAST == ICAST_C ) 
      { TABLEINFO_NPZ.PTRVAR[ITAB][IVAR] = (void*)((char*)PTRVAR) ; }

    TABLEINFO_NPZ.NVAR[ITAB]++ ;

  } // end of ivar

  return;

} // end  SNTABLE_ADDCOL_NPZ

// ==============================================
void SNTABLE_FILL_NPZ(int IDTABLE) {

  // Analog of SNTABLE_FILL_TEXT.
  // Append one row to the npz table for this IDTABLE: for each column,
  // copy the current value (read from the stored pointer) into the
  // column buffer, then mark the row complete. Compression and disk
  // spill happen inside npz_table_endrow every nflush rows.
  // Header is not needed here; npz_table_close writes the npy shapes.

  int  ITAB, NFILL, NVAR, IVAR, ICAST, HANDLE;
  char *VARNAME;
  char fnam[] = "SNTABLE_FILL_NPZ" ;

  // ----------- BEGIN ------------

  // find sparse index ITAB for this IDTABLE; abort if no table
  ITAB   = get_ITABLE_SPARSE(IDTABLE, fnam, 1); 
  NFILL  = TABLEINFO_NPZ.NFILL[ITAB];
  NVAR   = TABLEINFO_NPZ.NVAR[ITAB];
  HANDLE = TABLEINFO_NPZ.HANDLE[ITAB];

  if ( NVAR <= 0 ) {
    sprintf(MSGERR1, "No columns defined for IDTABLE=%d (%s)", 
	    IDTABLE, TABLEINFO_NPZ.TBNAME[ITAB]);
    sprintf(MSGERR2, "Call SNTABLE_ADDCOL_NPZ before SNTABLE_FILL_NPZ");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  // exactly one append per column per row, then one endrow
  for ( IVAR=0; IVAR < NVAR ; IVAR++ ) {
    ICAST   = TABLEINFO_NPZ.ICAST[ITAB][IVAR] ;
    VARNAME = TABLEINFO_NPZ.VARNAME[ITAB][IVAR] ;

    if ( ICAST < 0 || TABLEINFO_NPZ.PTRVAR[ITAB][IVAR] == NULL ) {
      sprintf(MSGERR1, "Undefined ICAST/pointer for IVAR=%d (%s)  IDTABLE=%d", 
	      IVAR, VARNAME, IDTABLE);
      sprintf(MSGERR2, "Something is messed up."); 
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
    }

    npz_table_append(HANDLE, TABLEINFO_NPZ.ICOL[ITAB][IVAR], 
		     TABLEINFO_NPZ.PTRVAR[ITAB][IVAR] );
  }

  npz_table_endrow(HANDLE);

  // increment number of FILL calls
  TABLEINFO_NPZ.NFILL[ITAB]++ ;  (void)NFILL;

} // end SNTABLE_FILL_NPZ

// ==============================================
static void npz_json_escape(std::string &out, const char *s) {
  // append s to out as a JSON string body (no surrounding quotes)
  for ( ; *s; s++ ) {
    unsigned char c = (unsigned char)*s;
    if      ( c == '"'  ) { out += "\\\""; }
    else if ( c == '\\' ) { out += "\\\\"; }
    else if ( c == '\n' ) { out += "\\n"; }
    else if ( c == '\t' ) { out += "\\t"; }
    else if ( c == '\r' ) { out += "\\r"; }
    else if ( c < 0x20  ) { char b[8]; sprintf(b,"\\u%04x",c); out += b; }
    else                  { out += (char)c; }       // UTF-8 bytes pass through
  }
}

std::string get_npz_header_json(int IDTABLE) {

  // Oct 2026
  // Return SNTABLE_CREATE_INFO.HEADER_DICT for this IDTABLE as a JSON object string.
  // Same information that SNTABLE_WRITE_HEADER_TEXT writes to the text header
  // (empty keys, which text writes as a blank '#' line, are skipped).
  // Key order is preserved. Values are always strings.
  //
  // read json header with python command:
  //   hdr = json.loads(data.zip.comment.decode())
  //

  char fnam[] = "get_npz_header_json" ;
  int  ITAB_SPARSE = get_ITABLE_SPARSE(IDTABLE, fnam, 1);   // index into SNTABLE_CREATE_INFO
  TABLE_HEADER_DICT_DEF *HEADER_DICT = &SNTABLE_CREATE_INFO.HEADER_DICT[ITAB_SPARSE];
  std::string json = "{";
  int j, nout = 0;

  for(j=0; j < HEADER_DICT->N_ITEM; j++ ) {
    const char *key = HEADER_DICT->KEY_STRING[j];
    const char *val = HEADER_DICT->VAL_STRING[j];
    if ( !key || strlen(key) == 0 ) { continue; }
    if ( nout > 0 ) { json += ", "; }
    json += "\""; npz_json_escape(json, key); json += "\": \"";
    npz_json_escape(json, val ? val : ""); json += "\"";
    nout++ ;
  }
  json += "}";
  return json;

} // end get_npz_header_json

void CLOSE_NPZFILE(void) {

  // Oct 2026: before closing each table, pass its header dictionary to the writer
  //           (stored as the zip archive comment).

  int ITAB;

  for(ITAB=0; ITAB < TABLEINFO_NPZ.NTABLE; ITAB++ ) {
    int HANDLE = TABLEINFO_NPZ.HANDLE[ITAB];
    std::string json = get_npz_header_json(TABLEINFO_NPZ.IDTABLE[ITAB]);
    npz_table_set_header(HANDLE, json.c_str());

    printf("   Close NPZ table file: %s \n", TABLEINFO_NPZ.FILENAME[ITAB]);
    fflush(stdout);
    npz_table_close(HANDLE);        // was npz_table_close(ITAB): wrong if HANDLE != ITAB
  }

  // reset so a second CLOSE call cannot re-close (and truncate) finished files
  TABLEINFO_NPZ.NTABLE = 0 ;

  return;

} // end CLOSE_NPZFILE
// ==============================================
int get_npz_filetype(char *npz_file) {

  // Created Sep 29 2026 by R.Kessler and Claude.
  // For npz file written by create_covariance.py,
  // return ISFILE_NPZ_COV *npz_file contains a 'cov' data type;
  // return ISFILE_NPZ_FACTORIEZE if *npz_file contains 'diag' data type.

  int ISFILE_TYPE = -9;
  char fnam[] = "get_npz_filetype" ;  (void)fnam;

  // ---------- BEGIN -----------

  cnpy::npz_t my_npz = cnpy::npz_load(npz_file);

  if ( my_npz.count("cov") > 0 ) {
    ISFILE_TYPE = ISFILE_NPZ_COV;
  }
  else if ( my_npz.count("diag") > 0 ) {
    ISFILE_TYPE = ISFILE_NPZ_FACTORIZE;
  }

  return ISFILE_TYPE;


} // end get_npz_filetype

int read_npz_covmat(char *npz_file, double *array1d) {

  // Created Apr 14 2024
  // Read cov matrix from npz file of the form
  //    nsn = nsn
  //    cov = cov(utria).float
  //
  // The cov in npz file is float (4 byte) and upper triangular half
  // to save disk space. 
  // However, the returned cov is the full N x N in double precision.

  int k0, k1, j1du=0, j1d_tmp, NSN;
  double cov;
  bool CHECK_COV_LOADING = false;   int *n_cov_load;
  char fnam[] = "read_npz_covmat" ;

  // ------------ BEGIN --------------

  // load the entire npz file
  cnpy::npz_t  my_npz = cnpy::npz_load(npz_file);

  // pick off NSN
  cnpy::NpyArray nsn_npz = my_npz["nsn"];
  NSN = *nsn_npz.data<int>();

  // read upper-triangle part of cov
  cnpy::NpyArray cov_utri_npz = my_npz["cov"];
  float *cov_utri = cov_utri_npz.data<float>() ;

  // - - - - - - -
  // load full NSN x NSN output array and convert float to double

  if ( CHECK_COV_LOADING ) { 
    n_cov_load = (int*) malloc( NSN*NSN * sizeof(int) );
    for(j1d_tmp = 0; j1d_tmp < NSN*NSN; j1d_tmp++ ) { n_cov_load[j1d_tmp] = false; }
  }

  for (k0=0; k0 < NSN; k0++ ) {
    for (k1=k0; k1 < NSN; k1++ ) {
      cov    = cov_utri[j1du];
       
      j1d_tmp = k0*NSN + k1;
      array1d[j1d_tmp] = (double)cov ;
      if ( CHECK_COV_LOADING ) { n_cov_load[j1d_tmp]++ ; }

      if ( k0 != k1 ) {
	j1d_tmp = k1*NSN + k0;
	array1d[j1d_tmp] = (double)cov ;
	if ( CHECK_COV_LOADING ) { n_cov_load[j1d_tmp]++ ; }
      }

      j1du+=1 ;
    }
  }


  // - - - - - - - -

  if ( CHECK_COV_LOADING ) {    
    for(j1d_tmp = 0; j1d_tmp < NSN*NSN; j1d_tmp++ ) {
      int n = n_cov_load[j1d_tmp] ;
      if ( n != 1 ) {
	printf(" xxx %s ERROR: j1d=%d  n_cov_load = %d \n",
	       fnam, j1d_tmp, n ); fflush(stdout);
      }
    }
  }

  // - - - -
  return NSN;

} // read_npz_array


// ===========================================================
int read_npz_factorized(char *npz_file, double **diag, double **U, int *K) {

  // Created Sep 2026
  // Read the additive factorized covariance product D + U U^T,
  // written by create_covariance.py's write_covariance_factorized().
  // npz file contains:
  //    nsn = NSN (int)
  //    diag = (NSN,) float64        statistical variance (diagonal of D)
  //    U    = (NSN,K) float64       per-systematic Delta_mu vectors
  // Both diag and U are malloc'ed here (caller must free); K is returned
  // via output arg (K=0 is valid: no systematics, e.g. NOSYS covopt).
  //
  // Returns NSN.

  char fnam[] = "read_npz_factorized" ;

  // ------------ BEGIN --------------

  cnpy::npz_t my_npz = cnpy::npz_load(npz_file);

  cnpy::NpyArray nsn_npz = my_npz["nsn"];
  int NSN = *nsn_npz.data<int>();

  cnpy::NpyArray diag_npz = my_npz["diag"];
  double *diag_data = diag_npz.data<double>();

  cnpy::NpyArray U_npz = my_npz["U"];
  int K_local = ( U_npz.shape.size() == 2 ) ? (int)U_npz.shape[1] : 0 ;
  double *U_data = U_npz.data<double>();

  if ( (int)diag_npz.shape[0] != NSN ) {
    // caller (wfit.c) does the fatal-abort sanity check on the return
    // value vs. expected NSN, matching the read_npz_covmat convention;
    // just warn here since sntools.h errmsg()/c1err are not available
    // in this translation unit.
    printf(" %s WARNING: diag array has %d elements, but nsn=%d (%s)\n",
	   fnam, (int)diag_npz.shape[0], NSN, npz_file); fflush(stdout);
  }

  *diag = (double*) malloc( NSN * sizeof(double) );
  memcpy(*diag, diag_data, NSN * sizeof(double) );

  long long n_u = (long long)NSN * (long long)K_local ;
  if ( n_u > 0 ) {
    *U = (double*) malloc( n_u * sizeof(double) );
    memcpy(*U, U_data, n_u * sizeof(double) );
  }
  else {
    *U = NULL ;
  }

  *K = K_local ;

  return NSN;

} // end read_npz_factorized



// =====================================================================
//   NPZ table writer: buffer rows, deflate + spill to temp files every
//   NFLUSH rows, assemble compressed .npz at close.
//
//   Handles tables > 4 GB (zip64 used only where needed).
//   Env overrides:  SNANA_NPZ_NFLUSH (default 100)
//                   SNANA_NPZ_LEVEL  (zlib 0-9, default 6; 0 = no compression)
// =====================================================================
#include <vector>
#include <string>
#include <stdexcept>
#include <stdint.h>
#include <zlib.h>

namespace {

  const int NPZ_ICAST_L=16, NPZ_ICAST_D=8, NPZ_ICAST_F=4,
            NPZ_ICAST_S=3,  NPZ_ICAST_I=2, NPZ_ICAST_C=1;   // match sntools_output.h

  struct NPZCol {
    std::string name;  int icast;  int wsize;
    std::vector<char> buf;                 // at most NFLUSH rows
    std::string tmpfile;  FILE *fp;        // compressed spill file
    z_stream *zs;                          // persistent deflate stream (lazy init)
    uLong     crc;                         // crc32 of uncompressed data
    uint64_t  nraw, ncomp;                 // uncompressed / compressed byte counts
    NPZCol() : icast(0), wsize(0), fp(NULL), zs(NULL), crc(0), nraw(0), ncomp(0) {}
  };

  struct NPZTab {
    std::string file;  std::vector<NPZCol> col;
    size_t nrow, nbuf, nflush;  int level;  bool used;
    std::string header_json;  bool has_header;   // table header dictionary, stored as JSON
    NPZTab() : nrow(0), nbuf(0), nflush(100), level(6), used(false), has_header(false) {}
  };
  NPZTab NPZ_TAB[MXTABLE_NPZ];

  std::vector<unsigned char> ZOUT(1<<16);  // shared deflate output buffer

  void put16(std::vector<char> &v, uint16_t x) { for(int i=0;i<2;i++) v.push_back((char)((x>>(8*i))&0xFF)); }
  void put32(std::vector<char> &v, uint32_t x) { for(int i=0;i<4;i++) v.push_back((char)((x>>(8*i))&0xFF)); }
  void putstr(std::vector<char> &v, const std::string &s) { v.insert(v.end(), s.begin(), s.end()); }

  // deflate n bytes into this column's spill file. flushflag = Z_NO_FLUSH or Z_FINISH.
  // Input is fed in slices <= 1 GB because zlib's avail_in is 32-bit.
  void col_deflate(NPZCol &c, int level, const char *in, size_t n, int flushflag) {
    if ( !c.zs ) {
      c.zs = new z_stream;  memset(c.zs, 0, sizeof(z_stream));
      // raw deflate (windowBits<0): zip stores raw streams.  memLevel 8 => ~256 KB/column
      if ( deflateInit2(c.zs, level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK )
        { throw std::runtime_error("npz: deflateInit2 failed"); }
      c.fp = fopen(c.tmpfile.c_str(), "wb");
      if ( !c.fp ) { throw std::runtime_error("npz: cannot open " + c.tmpfile); }
      c.crc = crc32(0L, Z_NULL, 0);
    }
    const size_t SLICE = (size_t)1 << 30;
    size_t done = 0;
    do {
      size_t m = n - done;  if ( m > SLICE ) { m = SLICE; }
      bool last = ( done + m >= n );
      if ( m > 0 ) { c.crc = crc32(c.crc, (const Bytef*)(in+done), (uInt)m); }
      c.nraw += m;

      c.zs->next_in  = (Bytef*)(in + done);
      c.zs->avail_in = (uInt)m;
      int ff = last ? flushflag : Z_NO_FLUSH;
      do {
        c.zs->next_out  = ZOUT.data();
        c.zs->avail_out = (uInt)ZOUT.size();
        int ret = deflate(c.zs, ff);
        if ( ret == Z_STREAM_ERROR ) { throw std::runtime_error("npz: deflate error"); }
        size_t have = ZOUT.size() - c.zs->avail_out;
        if ( have && fwrite(ZOUT.data(), 1, have, c.fp) != have )
          { throw std::runtime_error("npz: write error " + c.tmpfile); }
        c.ncomp += have;
      } while ( c.zs->avail_out == 0 );
      done += m;
    } while ( done < n );
  }

  void npz_flush(NPZTab &t) {
    if ( t.nbuf == 0 ) { return; }
    for (size_t i=0; i < t.col.size(); i++) {
      col_deflate(t.col[i], t.level, t.col[i].buf.data(), t.col[i].buf.size(), Z_NO_FLUSH);
      t.col[i].buf.clear();
    }
    t.nbuf = 0;
  }

  // Build .npy v1.0 header for a 1D fixed-width byte-string array, dtype '|S<wsize>'.
  // cnpy has no string dtype (char maps to int8), so build it by hand.
  std::vector<char> npy_header_str(size_t nrow, int wsize) {
    std::string d = "{'descr': '|S" + std::to_string(wsize) +
                    "', 'fortran_order': False, 'shape': (" +
                    std::to_string(nrow) + ",), }";
    size_t total = 10 + d.size() + 1;           // magic(6)+version(2)+len(2) + dict + '\\n'
    d.append((64 - total % 64) % 64, ' ');      // pad so data starts 64-byte aligned
    d.push_back('\n');
    std::vector<char> h = { '\x93','N','U','M','P','Y', 1, 0 };
    put16(h, (uint16_t)d.size());
    putstr(h, d);
    return h;
  }

  std::vector<char> npy_header_for(int icast, size_t nrow, int wsize) {
    std::vector<size_t> s1 = { nrow };
    switch (icast) {
      case NPZ_ICAST_D: return cnpy::create_npy_header<double>(s1);
      case NPZ_ICAST_F: return cnpy::create_npy_header<float>(s1);
      case NPZ_ICAST_I: return cnpy::create_npy_header<int>(s1);
      case NPZ_ICAST_S: return cnpy::create_npy_header<short>(s1);
      case NPZ_ICAST_L: return cnpy::create_npy_header<long long>(s1);
      default:          return npy_header_str(nrow, wsize);         // fixed-width strings
    }
  }

  void put64(std::vector<char> &v, uint64_t x) { for(int i=0;i<8;i++) v.push_back((char)((x>>(8*i))&0xFF)); }
} // end anonymous namespace


int npz_table_create(char *fileName) {

  int LDMP = 0;
  char fnam[] = "npz_table_create";

  for (int h=0; h < MXTABLE_NPZ; h++) {
    if ( !NPZ_TAB[h].used ) {
      NPZ_TAB[h] = NPZTab();  NPZ_TAB[h].used = true;  NPZ_TAB[h].file = fileName;
      const char *e;
      if ( (e = getenv("SNANA_NPZ_NFLUSH")) && atoi(e) > 0 ) 
	{ NPZ_TAB[h].nflush = atoi(e); }
      if ( (e = getenv("SNANA_NPZ_LEVEL"))  && atoi(e) >= 0 && atoi(e) <= 9 ) 
	{ NPZ_TAB[h].level = atoi(e); }

      if(LDMP) { printf(" xxx %s: create %s for h=%d \n", fnam, fileName, h); fflush(stdout); }

      return h;
    }
  }
  return -1;
}

void npz_table_set_header(int h, const char *json) {
  // Store header dictionary (as a JSON string) for table h.
  // Zip comment length is a 16-bit field, so the JSON is limited to 65535 bytes.
  // It is written as the zip archive comment (not a table column) when the table is closed.
  if ( h < 0 || h >= MXTABLE_NPZ || !NPZ_TAB[h].used || !json ) { return; }
  NPZ_TAB[h].header_json = json;
  NPZ_TAB[h].has_header  = ( NPZ_TAB[h].header_json.size() > 0 );
  if ( NPZ_TAB[h].header_json.size() > 0xFFFF ) {
    printf(" WARNING npz_table_set_header: header JSON is %lu bytes (max 65535); not written to %s\n",
	   (unsigned long)NPZ_TAB[h].header_json.size(), NPZ_TAB[h].file.c_str());
    NPZ_TAB[h].has_header = false;
  }
}

int npz_table_addcol(int h, char *varName, int icast, int strlen_max) {

  // ======================================================================
  // npz_table_addcol
  //
  // FREQUENCY: called ONCE PER COLUMN, during table setup (after
  //            npz_table_create, before the first row is written).
  //            It is not part of the per-event loop.
  //
  // Defines one column of table 'h' and returns its index (0,1,2,...),
  // which the caller must save and pass as 'icol' to npz_table_append.
  //
  // Inputs:
  //   h          : table handle returned by npz_table_create
  //   varName    : column name; becomes "<varName>.npy" inside the .npz
  //   icast      : data type, one of NPZ_ICAST_[D,F,I,S,L,C]
  //                (double, float, int, short, long long, char-string)
  //   strlen_max : fixed string width in bytes; used ONLY for NPZ_ICAST_C,
  //                ignored for numeric types
  //
  // Returns: column index (>= 0), or -1 if icast is not recognized
  //          (in which case no column is added).
  //
  // Notes:
  //  - Only metadata is created here. No file is opened and no zlib stream
  //    is allocated; both are created lazily in col_deflate() on the first
  //    flush. Memory cost is just an empty NPZCol until rows arrive.
  //  - Add ALL columns before the first npz_table_endrow. A column added
  //    after rows have been flushed would be shorter than the others, and
  //    the npz file would be inconsistent. There is no check for this.
  // ======================================================================

  char fnam[] = "npz_table_addcol";  (void)fnam;

  NPZCol c;  c.name = varName;  c.icast = icast;
  switch (icast) {
    case NPZ_ICAST_D: c.wsize = sizeof(double);        break;
    case NPZ_ICAST_F: c.wsize = sizeof(float);         break;
    case NPZ_ICAST_I: c.wsize = sizeof(int);           break;
    case NPZ_ICAST_S: c.wsize = sizeof(short int);     break;
    case NPZ_ICAST_L: c.wsize = sizeof(long long int); break;
    case NPZ_ICAST_C: c.wsize = strlen_max;            break;
    default: return -1;
  }
  c.tmpfile = NPZ_TAB[h].file + ".tmp." + std::to_string(NPZ_TAB[h].col.size()) + ".z";
  NPZ_TAB[h].col.reserve(1024);   // never reallocate after streams exist
  NPZ_TAB[h].col.push_back(c);
  return (int)NPZ_TAB[h].col.size() - 1;
}

void npz_table_append(int h, int icol, void *ptr) {

  // ======================================================================
  // npz_table_append
  //
  // FREQUENCY: called ONCE PER COLUMN PER ROW (i.e. per event). For each
  //            event, call it once for every column, then call
  //            npz_table_endrow once. Total calls = NCOL * NROW.
  //
  // Copies the current value of one column into that column's in-memory
  // buffer. The data are compressed and written to disk later, by
  // npz_flush, which npz_table_endrow triggers every nflush rows.
  //
  // Inputs:
  //   h    : table handle
  //   icol : column index returned by npz_table_addcol
  //   ptr  : address of the value for THIS row. It is read immediately,
  //          so the pointer need not stay valid after the call. It must
  //          point to a type matching the column's icast:
  //            D -> double*, F -> float*, I -> int*, S -> short*,
  //            L -> long long*, C -> char* (null-terminated string)
  //
  // Notes:
  //  - Numeric types: wsize bytes are copied verbatim.
  //  - Strings: the slot is zero-filled first, then at most wsize bytes
  //    are copied (stopping at the string's NUL). Strings are therefore
  //    NUL-padded to strlen_max. A string of length >= strlen_max is
  //    truncated and has no terminator, as in numpy fixed-width 'S' arrays.
  //  - Nothing is checked, so every column must receive exactly one append
  //    per row. Missing or extra appends silently desynchronize column lengths.
  //  - Memory use is bounded: a column buffer holds at most nflush rows
  //    (default 100) before it is cleared by npz_flush.
  // =========================================================

  NPZCol &c = NPZ_TAB[h].col[icol];
  size_t n0 = c.buf.size();
  c.buf.resize(n0 + c.wsize, 0);
  if ( c.icast == NPZ_ICAST_C ) { memcpy(&c.buf[n0], ptr, strnlen((char*)ptr, c.wsize)); }
  else                          { memcpy(&c.buf[n0], ptr, c.wsize); }
}

void npz_table_putcol(int h, int icol, void *ptr, size_t nrow) {

  // ======================================================================
  // npz_table_putcol   (Nov 2026)
  //
  // FREQUENCY: called ONCE PER COLUMN (not per row) to write an entire column 
  //            that is already in memory. Use this INSTEAD of 
  //            npz_table_append/endrow for a table; do not mix the two methods 
  //            in one table.
  //
  //   h    : table handle
  //   icol : column index returned by npz_table_addcol
  //   ptr  : array of nrow elements matching the column's icast
  //          (strings: flat array of nrow*strlen_max bytes, NUL-padded)
  //   nrow : number of rows; must be the same for every column.
  //
  // The data are compressed directly to the column's temp file; no 
  // copy of the column is made. All columns must be added (npz_table_addcol) 
  // before the first putcol.
  // ======================================================================

  NPZTab &t = NPZ_TAB[h];
  NPZCol &c = t.col[icol];

  if ( t.nbuf > 0 || (c.nraw > 0) ) 
    { throw std::runtime_error("npz: putcol called twice for column " + c.name + " or after row-wise append"); }
  if ( icol > 0 && t.nrow != nrow ) 
    { throw std::runtime_error("npz: putcol nrow mismatch for column " + c.name); }
  t.nrow = nrow ;

  // when nrow=0, still create the (empty) stream so that the column is closed properly
  col_deflate(c, t.level, (const char*)ptr, nrow * (size_t)c.wsize, Z_NO_FLUSH);
}

void npz_table_endrow(int h) {

  // ======================================================================
  // npz_table_endrow
  //
  // FREQUENCY: called ONCE PER ROW (i.e. per event), after
  //            npz_table_append has been called for every column of that
  //            row. It is not per column and not per table.
  //
  // Marks the current row as complete and decides whether to spill to disk:
  //   - nrow : running total of rows. It is used at close time to write
  //            the .npy shape header (nrow) for every column.
  //   - nbuf : rows currently held in memory since the last flush.
  //   - every nflush rows (default 100; env SNANA_NPZ_NFLUSH overrides),
  //     npz_flush() deflates each column buffer into its temp file and
  //     empties the buffers, so memory stays bounded however many events
  //     are written.
  //
  // Any rows still buffered when the table ends (fewer than nflush) are
  // written by npz_table_close, which finishes each deflate stream.
  //
  // Per-event pattern:
  //   for each event {
  //     for each column icol:  npz_table_append(h, icol, &value[icol]);
  //     npz_table_endrow(h);
  //   }
  // ===================================================
  NPZTab &t = NPZ_TAB[h];
  t.nrow++;  t.nbuf++;
  if ( t.nbuf >= t.nflush ) { npz_flush(t); }
}

void npz_table_close(int h) {

  // Writes a zip archive.  Entries / end-of-file records use the zip64 format
  // only when a 32-bit field would overflow, so small tables stay plain zip32
  // (readable by any unzip) while tables > 4 GB work too (numpy/zipfile read zip64).
  // On any failure, temp files and the partial .npz are removed, then the
  // exception propagates.

  // ---------- BEGIN ----------

  NPZTab &t = NPZ_TAB[h];
  FILE *fz = NULL;
  // threshold at which zip64 records are used (compile with -DNPZ_ZIP64_LIMIT=<n> to test
  // zip64 paths on small tables; 0xFFFFFFFF markers in 32-bit fields are always correct)
#ifndef NPZ_ZIP64_LIMIT
#define NPZ_ZIP64_LIMIT 0xFFFFFFFFULL
#endif
  const uint64_t M32 = NPZ_ZIP64_LIMIT;

  try {
    fz = fopen(t.file.c_str(), "wb");
    if ( !fz ) { throw std::runtime_error("npz: cannot open " + t.file); }

    std::vector<char> central;
    uint64_t offset = 0;                     // running offset in output file
    std::vector<char> cp(1<<20);

    for (size_t i=0; i < t.col.size(); i++) {
      NPZCol &c = t.col[i];

      // finish deflate stream with whatever is left in the buffer
      col_deflate(c, t.level, c.buf.data(), c.buf.size(), Z_FINISH);
      fclose(c.fp);  c.fp = NULL;
      deflateEnd(c.zs);  delete c.zs;  c.zs = NULL;
      c.buf.clear();  c.buf.shrink_to_fit();

      // npy header is only known now (needs nrow). Put it first in the deflate stream as a
      // non-final STORED block (0x00, LEN, ~LEN, bytes); the column's own stream follows.
      std::vector<char> hdr = npy_header_for(c.icast, t.nrow, c.wsize);
      uint16_t hl = (uint16_t)hdr.size();
      uLong crc_hdr = crc32(crc32(0L, Z_NULL, 0), (const Bytef*)hdr.data(), (uInt)hdr.size());
      uint32_t crc  = (uint32_t)crc32_combine(crc_hdr, c.crc, (z_off_t)c.nraw);
      uint64_t usize = hdr.size() + c.nraw;
      uint64_t csize = 5 + hdr.size() + c.ncomp;

      bool big_sizes  = ( usize >= M32 || csize >= M32 );   // needs zip64 extra in local+central
      bool big_offset = ( offset >= M32 );                  // needs zip64 extra in central only

      std::string fname = c.name + ".npy";

      // ---- local header ----
      std::vector<char> lh;
      put32(lh, 0x04034b50);
      put16(lh, big_sizes ? 45 : 20);  put16(lh, 0);  put16(lh, 8);      // method 8 = deflate
      put16(lh, 0); put16(lh, 0);  put32(lh, crc);
      put32(lh, big_sizes ? 0xFFFFFFFFu : (uint32_t)csize);
      put32(lh, big_sizes ? 0xFFFFFFFFu : (uint32_t)usize);
      put16(lh, (uint16_t)fname.size());
      put16(lh, big_sizes ? 20 : 0);                                      // extra-field length
      putstr(lh, fname);
      if ( big_sizes ) {                  // zip64 extra: local header needs BOTH sizes (usize first)
        put16(lh, 0x0001);  put16(lh, 16);  put64(lh, usize);  put64(lh, csize);
      }

      if ( fwrite(lh.data(), 1, lh.size(), fz) != lh.size() )
        { throw std::runtime_error("npz: write error " + t.file); }
      char sb[5] = { 0, (char)(hl & 0xFF), (char)(hl >> 8), (char)(~hl & 0xFF), (char)((~hl >> 8) & 0xFF) };
      fwrite(sb, 1, 5, fz);
      fwrite(hdr.data(), 1, hdr.size(), fz);

      FILE *fi = fopen(c.tmpfile.c_str(), "rb");
      if ( !fi ) { throw std::runtime_error("npz: cannot reopen " + c.tmpfile); }
      size_t n;  uint64_t ncopy = 0;
      while ( (n = fread(cp.data(), 1, cp.size(), fi)) > 0 ) {
        if ( fwrite(cp.data(), 1, n, fz) != n ) { fclose(fi); throw std::runtime_error("npz: write error " + t.file); }
        ncopy += n;
      }
      fclose(fi);
      if ( ncopy != c.ncomp ) { throw std::runtime_error("npz: temp file size mismatch " + c.tmpfile); }
      remove(c.tmpfile.c_str());

      // ---- central-directory entry ----
      // zip64 extra holds only the fields whose 32-bit slot is 0xFFFFFFFF, in this order:
      // usize, csize, local-header offset.
      std::vector<char> extra;
      if ( big_sizes )  { put64(extra, usize);  put64(extra, csize); }
      if ( big_offset ) { put64(extra, offset); }
      bool z64 = ( big_sizes || big_offset );

      put32(central, 0x02014b50);
      put16(central, z64 ? 45 : 20);  put16(central, z64 ? 45 : 20);  put16(central, 0);
      put16(central, 8); put16(central, 0); put16(central, 0); put32(central, crc);
      put32(central, big_sizes  ? 0xFFFFFFFFu : (uint32_t)csize);
      put32(central, big_sizes  ? 0xFFFFFFFFu : (uint32_t)usize);
      put16(central, (uint16_t)fname.size());
      put16(central, z64 ? (uint16_t)(4 + extra.size()) : 0);
      put16(central, 0);  put16(central, 0);  put16(central, 0);  put32(central, 0);
      put32(central, big_offset ? 0xFFFFFFFFu : (uint32_t)offset);
      putstr(central, fname);
      if ( z64 ) { put16(central, 0x0001);  put16(central, (uint16_t)extra.size());
                   central.insert(central.end(), extra.begin(), extra.end()); }

      offset += lh.size() + 5 + hdr.size() + c.ncomp;
    }

    // ---- end of central directory (zip64 records only if needed) ----
    uint64_t cd_offset = offset,  cd_size = central.size(),  ncol = t.col.size();
    bool eocd64 = ( cd_offset >= M32 || cd_size >= M32 || ncol >= 0xFFFF );

    std::vector<char> tail;
    if ( eocd64 ) {
      put32(tail, 0x06064b50);  put64(tail, 44);              // zip64 EOCD record
      put16(tail, 45);  put16(tail, 45);
      put32(tail, 0);   put32(tail, 0);
      put64(tail, ncol);  put64(tail, ncol);
      put64(tail, cd_size);  put64(tail, cd_offset);
      put32(tail, 0x07064b50);  put32(tail, 0);               // zip64 EOCD locator
      put64(tail, cd_offset + cd_size);  put32(tail, 1);
    }
    put32(tail, 0x06054b50);  put16(tail, 0);  put16(tail, 0);
    put16(tail, ncol >= 0xFFFF ? 0xFFFF : (uint16_t)ncol);
    put16(tail, ncol >= 0xFFFF ? 0xFFFF : (uint16_t)ncol);
    put32(tail, cd_size   >= M32 ? 0xFFFFFFFFu : (uint32_t)cd_size);
    put32(tail, cd_offset >= M32 ? 0xFFFFFFFFu : (uint32_t)cd_offset);
    // archive comment: header dictionary as JSON. Not a zip entry, so it never appears in npz.files.
    // python: json.loads(np.load(f).zip.comment.decode())
    if ( t.has_header ) { put16(tail, (uint16_t)t.header_json.size());  putstr(tail, t.header_json); }
    else               { put16(tail, 0); }

    if ( fwrite(central.data(), 1, central.size(), fz) != central.size() ||
         fwrite(tail.data(),    1, tail.size(),    fz) != tail.size() )
      { throw std::runtime_error("npz: write error " + t.file); }
    if ( fclose(fz) != 0 ) { fz = NULL; throw std::runtime_error("npz: close error " + t.file); }
    fz = NULL;
  }
  catch (...) {
    // cleanup: close streams/files, delete temp files and partial output
    if ( fz ) { fclose(fz); }
    for (size_t i=0; i < t.col.size(); i++) {
      NPZCol &c = t.col[i];
      if ( c.zs ) { deflateEnd(c.zs); delete c.zs; c.zs = NULL; }
      if ( c.fp ) { fclose(c.fp); c.fp = NULL; }
      remove(c.tmpfile.c_str());
    }
    remove(t.file.c_str());
    t = NPZTab();
    throw;
  }

  t = NPZTab();
} // end npz_table_close


// ==========================================================================
//   Nov 26 2026: NPZ implementation of SNTABLE_APPENDCOL_EXEC
//   (generic part is in sntools_output.cpp)
//
//   Column-oriented: each input column is inflated once, the rows to keep
//   are gathered, and the whole output column is handed to the npz writer
//   (npz_table_putcol). There is no row-by-row I/O. Rows are matched to the
//   stored columns with a hash table, so processing time is linear in NROW.
//   Peak memory is about one input column plus one output column.
// ==========================================================================

namespace {

  // element access for an input npz column
  inline bool npz_is_str(const npzrd::COL &c) { return (c.kind=='S' || c.kind=='U'); }

  inline double npz_get_d(const npzrd::COL &c, const char *p) {
    if ( c.kind=='f' && c.wsize==8 ) { double x; memcpy(&x,p,8); return x; }
    if ( c.kind=='f' )               { float  x; memcpy(&x,p,4); return (double)x; }
    if ( c.wsize==8 ) { if (c.kind=='i') { long long x; memcpy(&x,p,8); return (double)x; } 
                        else { unsigned long long x; memcpy(&x,p,8); return (double)x; } }
    if ( c.wsize==4 ) { if (c.kind=='i') { int x; memcpy(&x,p,4); return (double)x; } 
                        else { unsigned x; memcpy(&x,p,4); return (double)x; } }
    if ( c.wsize==2 ) { if (c.kind=='i') { short x; memcpy(&x,p,2); return (double)x; } 
                        else { unsigned short x; memcpy(&x,p,2); return (double)x; } }
    return (c.kind=='i') ? (double)(signed char)p[0] : (double)(unsigned char)p[0];
  }

  inline long long npz_get_l(const npzrd::COL &c, const char *p) {
    if ( c.kind=='f' ) { return (long long) npz_get_d(c,p); }
    if ( c.wsize==8 ) { long long x; memcpy(&x,p,8); return x; }
    if ( c.wsize==4 ) { if (c.kind=='i') { int x; memcpy(&x,p,4); return x; } 
                        else { unsigned x; memcpy(&x,p,4); return x; } }
    if ( c.wsize==2 ) { if (c.kind=='i') { short x; memcpy(&x,p,2); return x; } 
                        else { unsigned short x; memcpy(&x,p,2); return x; } }
    return (c.kind=='i') ? (long long)(signed char)p[0] : (long long)(unsigned char)p[0];
  }

  // string element -> s (NUL terminated, trailing blanks removed, max n-1 chars)
  inline void npz_get_s(const npzrd::COL &c, const char *p, char *s, int n) {
    int k=0;
    if ( npz_is_str(c) ) {
      int nch = (c.kind=='S') ? c.wsize : c.wsize/4 ;
      for(int j=0; j < nch && k < n-1; j++ ) {
	char ch = (c.kind=='S') ? p[j] : p[4*j];
	if ( ch == 0 ) { break; }
	s[k++] = ch;
      }
      while ( k > 0 && (s[k-1]==' ' || s[k-1]=='\t') ) { k--; }
      s[k] = 0;
    }
    else if ( c.kind == 'f' ) { snprintf(s, n, "%.10g", npz_get_d(c,p)); }
    else                      { snprintf(s, n, "%lld",  npz_get_l(c,p)); }
  }

  // write one value into output array element at dst, with output cast ocast and width ow.
  // Source is either an input npz element (c,p) or a double/string (for fill values).
  inline void npz_put_num(char *dst, int ocast, double d, long long l, bool isint) {
    if      ( ocast == ICAST_D ) { double x = isint ? (double)l : d;  memcpy(dst, &x, 8); }
    else if ( ocast == ICAST_F ) { float  x = isint ? (float)l  : (float)d;  memcpy(dst, &x, 4); }
    else if ( ocast == ICAST_I ) { int    x = isint ? (int)l    : (int)d;    memcpy(dst, &x, 4); }
    else if ( ocast == ICAST_S ) { short  x = isint ? (short)l  : (short)d;  memcpy(dst, &x, 2); }
    else if ( ocast == ICAST_L ) { long long x = isint ? l      : (long long)d;  memcpy(dst, &x, 8); }
  }

  inline void npz_put_from_src(char *dst, int ocast, int ow, const npzrd::COL &c, const char *p) {
    if ( ocast == ICAST_C ) {
      char s[256];  npz_get_s(c, p, s, (int)sizeof(s));
      int n = (int)strlen(s);  if ( n > ow ) { n = ow; }
      memset(dst, 0, ow);  memcpy(dst, s, n);
    }
    else if ( npz_is_str(c) ) {
      char s[256];  npz_get_s(c, p, s, (int)sizeof(s));
      npz_put_num(dst, ocast, strtod(s,NULL), strtoll(s,NULL,10), false);
    }
    else if ( c.kind == 'f' ) { npz_put_num(dst, ocast, npz_get_d(c,p), 0, false); }
    else                      { npz_put_num(dst, ocast, 0.0, npz_get_l(c,p), true); }
  }

  // fill value (-9 for column missing from a file, -888 for undefined row)
  inline void npz_put_fill(char *dst, int ocast, int ow, int val) {
    if ( ocast == ICAST_C ) { 
      char s[16];  sprintf(s, "%d", val);
      int n = (int)strlen(s);  if ( n > ow ) { n = ow; }
      memset(dst, 0, ow);  memcpy(dst, s, n);
    }
    else { npz_put_num(dst, ocast, (double)val, (long long)val, true); }
  }

  // element from stored column icol, row irow
  inline void npz_put_from_stored(char *dst, int ocast, int ow, int icol, int irow) {
    SNTABLE_APPENDCOL_DEF *A = &SNTABLE_APPENDCOL ;
    int ICAST = A->ICAST[icol] ;
    void *P   = A->PTR[icol] ;
    if      ( ICAST == ICAST_D ) { memcpy(dst, (double*)P + irow, 8); }
    else if ( ICAST == ICAST_F ) { memcpy(dst, (float*)P + irow, 4); }
    else if ( ICAST == ICAST_I ) { memcpy(dst, (int*)P + irow, 4); }
    else if ( ICAST == ICAST_S ) { memcpy(dst, (short int*)P + irow, 2); }
    else if ( ICAST == ICAST_L ) { memcpy(dst, (long long int*)P + irow, 8); }
    else if ( ICAST == ICAST_C ) {
      int isize = A->ISIZE[icol];
      char *c   = (char*)P + (size_t)irow*isize ;
      int  n    = (int)strnlen(c, isize);  if ( n > ow ) { n = ow; }
      memset(dst, 0, ow);  memcpy(dst, c, n);
    }
    (void)ocast;
  }

  inline int npz_wsize_cast(int icast, int strlen_max) {
    switch(icast) {
    case ICAST_D: return 8;  case ICAST_F: return 4;  case ICAST_I: return 4;
    case ICAST_S: return 2;  case ICAST_L: return 8;
    default:      return (strlen_max > 0) ? strlen_max : 1 ;
    }
  }

  // promote cast when same column has different types in different input files
  int npz_promote_cast(int a, int b) {
    if ( a < 0 ) { return b; }
    if ( a == b ) { return a; }
    if ( a == ICAST_C || b == ICAST_C ) { return ICAST_C; }
    if ( a == ICAST_D || b == ICAST_D ) { return ICAST_D; }
    bool aflt = (a == ICAST_F),  bflt = (b == ICAST_F);
    if ( aflt || bflt ) { return ICAST_D; }          // float mixed with integer
    if ( a == ICAST_L || b == ICAST_L ) { return ICAST_L; }
    if ( a == ICAST_I || b == ICAST_I ) { return ICAST_I; }
    return ICAST_S;
  }

  int npz_find_idcol(const npzrd::INFO &R, const char *VARNAME_ID) {
    const char *idlist[] = VARNAME_TABLE_IDLIST ;
    for(size_t i=0; i < R.col.size(); i++ ) 
      { if ( strcasecmp(R.col[i].name.c_str(), VARNAME_ID) == 0 ) { return (int)i; } }
    if ( strcmp(VARNAME_ID,"CID")==0 || strcmp(VARNAME_ID,"CCID")==0 ) {
      for(int j=0; j < N_VARNAME_TABLE_IDLIST; j++ ) 
	for(size_t i=0; i < R.col.size(); i++ ) 
	  { if ( strcasecmp(R.col[i].name.c_str(), idlist[j]) == 0 ) { return (int)i; } }
    }
    return -1;
  }

} // end anonymous namespace

// ----------------------------------------------------
void SNTABLE_APPENDCOL_EXEC_NPZ(void) {

  SNTABLE_APPENDCOL_DEF *A = &SNTABLE_APPENDCOL ;
  APPENDCOL_HASH_DEF  H ;
  int   NLIST, NFILE = A->NFILE_INP ;
  int   ifile, ivar, icol, NVAR_OUT, NCOL_OUT, NROW_INP=0, NMATCH=0, NUNMATCHED=0 ;
  bool  ALLROWS = ( A->OPTMASK == OPTMASK_APPENDCOL_ALLROWS );
  char  fnam[] = "SNTABLE_APPENDCOL_EXEC_NPZ" ;

  // -------------- BEGIN --------------

  NLIST = appendcol_check_stored(fnam);
  memset(&H, 0, sizeof(H));
  appendcol_hash_init(&H);

  // ---- open each input file (read metadata only) and keep info for all of them 
  std::vector<npzrd::INFO> F(NFILE);
  std::vector<FILE*> FP(NFILE, (FILE*)NULL);
  std::vector<int>  IVAR_ID(NFILE);
  std::string err;

  for(ifile=0; ifile < NFILE; ifile++ ) {
    if ( !npzrd::open_file(A->INPFILE[ifile], err) ) {
      print_preAbort_banner(fnam);
      sprintf(MSGERR1, "Could not read npz table file: %.150s", A->INPFILE[ifile]);
      sprintf(MSGERR2, "%.*s", 200, err.c_str());
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
    }
    F[ifile] = npzrd::R ;
    IVAR_ID[ifile] = npz_find_idcol(F[ifile], A->VARNAME_ID);
    if ( IVAR_ID[ifile] < 0 ) {
      snprintf(MSGERR1, sizeof(MSGERR1), "Could not find ID column '%.150s' in", A->VARNAME_ID);
      sprintf(MSGERR2, "%.*s", 200, A->INPFILE[ifile]);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
    }
    FP[ifile] = fopen(A->INPFILE[ifile], "rb");
    if ( !FP[ifile] ) {
      sprintf(MSGERR1, "Cannot open"); 
      sprintf(MSGERR2, "%.200s", A->INPFILE[ifile]);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
    }
    NROW_INP += (int)F[ifile].col[0].nrow ;
  }
  npzrd::R = npzrd::INFO();   // do not leave a table 'opened' for the table reader

  // ---- carried-over columns: user list, or all columns in first file
  std::vector<std::string> VARNAME_OUT;
  if ( A->NVAR_SELECT > 0 ) 
    { for(ivar=0; ivar < A->NVAR_SELECT; ivar++ ) { VARNAME_OUT.push_back(A->VARNAME_SELECT[ivar]); } }
  else
    { for(size_t i=0; i < F[0].col.size(); i++ ) { VARNAME_OUT.push_back(F[0].col[i].name); } }
  NVAR_OUT = (int)VARNAME_OUT.size();

  std::vector<char*> VARNAME_PTR(NVAR_OUT);
  for(ivar=0; ivar < NVAR_OUT; ivar++ ) { VARNAME_PTR[ivar] = (char*)VARNAME_OUT[ivar].c_str(); }

  std::vector<int> ICOL_REPLACE(NVAR_OUT), IS_APPENDED(A->NCOL+1);
  appendcol_replaced_index(NVAR_OUT, &VARNAME_PTR[0], &ICOL_REPLACE[0], &IS_APPENDED[0]);

  // IMAP[ifile][ivar] = column index in file, or -1
  std::vector< std::vector<int> > IMAP(NFILE, std::vector<int>(NVAR_OUT, -1));
  for(ifile=0; ifile < NFILE; ifile++ ) 
    for(ivar=0; ivar < NVAR_OUT; ivar++ ) 
      for(size_t i=0; i < F[ifile].col.size(); i++ ) 
	if ( F[ifile].col[i].name == VARNAME_OUT[ivar] ) { IMAP[ifile][ivar] = (int)i; break; }

  // ---- pass 1: match rows of each file to stored rows; decide rows to keep
  std::vector< std::vector<int> > ROWS(NFILE) ;   // input rows to keep, per file
  std::vector<int>  ILIST_OUT ;                   // stored index for each output row (or -1)
  std::vector<char> USED(NLIST+1, 0);
  std::vector<char> buf;
  
  for(ifile=0; ifile < NFILE; ifile++ ) {
    npzrd::COL &c = F[ifile].col[IVAR_ID[ifile]] ;
    size_t nrow = c.nrow ;
    if ( !npzrd::inflate_col(FP[ifile], c, buf, 0) || buf.size() < c.npy_off + (size_t)c.wsize*nrow ) {
      sprintf(MSGERR1, "Could not read ID column '%s'", c.name.c_str());
      sprintf(MSGERR2, "%.200s", A->INPFILE[ifile]);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
    }
    const char *base = &buf[0] + c.npy_off ;
    char ID[80];
    for(size_t r=0; r < nrow; r++ ) {
      npz_get_s(c, base + r*c.wsize, ID, (int)sizeof(ID));
      int m = appendcol_hash_find(&H, ID, ifile);
      if ( m < 0 && !ALLROWS ) { continue; }
      if ( m >= 0 ) { NMATCH++ ; USED[m] = 1; }
      ROWS[ifile].push_back((int)r);
      ILIST_OUT.push_back(m);
    }
  }
  size_t NROW_OUT = ILIST_OUT.size();
  for(int j=0; j < NLIST; j++ ) { if ( !USED[j] ) { NUNMATCHED++ ; } }

  // ---- define output columns: carried-over (maybe replaced), then appended
  struct OUTCOL { std::string name; int icast; int wsize; int src_stored; bool carried; int ivar; };
  std::vector<OUTCOL> OC;

  for(ivar=0; ivar < NVAR_OUT; ivar++ ) {
    OUTCOL o;  o.name = VARNAME_OUT[ivar];  o.carried = true;  o.ivar = ivar;
    o.src_stored = ICOL_REPLACE[ivar];
    if ( o.src_stored >= 0 ) {
      o.icast = A->ICAST[o.src_stored];
      o.wsize = npz_wsize_cast(o.icast, A->ISIZE[o.src_stored]);
    }
    else {
      int ic = -9, w = 1;
      for(ifile=0; ifile < NFILE; ifile++ ) {
	int i = IMAP[ifile][ivar];
	if ( i < 0 ) { continue; }
	npzrd::COL &c = F[ifile].col[i];
	ic = npz_promote_cast(ic, c.icast);
	int cw = (c.kind=='U') ? c.wsize/4 : c.wsize ;
	if ( npz_is_str(c) ) { if ( cw > w ) { w = cw; } } else { if ( 24 > w ) { w = 24; } }
      }
      if ( ic < 0 ) { ic = ICAST_D; }       // column is in no file; all values are -9
      o.icast = ic ;
      o.wsize = npz_wsize_cast(ic, w);
    }
    OC.push_back(o);
  }
  for(icol=0; icol < A->NCOL; icol++ ) {
    if ( !IS_APPENDED[icol] ) { continue; }
    OUTCOL o;  o.name = A->VARNAME[icol];  o.carried = false;  o.ivar = -1;
    o.src_stored = icol ;  o.icast = A->ICAST[icol];  o.wsize = npz_wsize_cast(o.icast, A->ISIZE[icol]);
    OC.push_back(o);
  }
  NCOL_OUT = (int)OC.size();

  // ---- create output table and define columns
  int h = npz_table_create(A->OUTFILE);
  if ( h < 0 ) {
    sprintf(MSGERR1, "Could not create npz table (too many open tables?)");
    sprintf(MSGERR2, "%.200s", A->OUTFILE);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }
  for(int k=0; k < NCOL_OUT; k++ ) {
    int ic = npz_table_addcol(h, (char*)OC[k].name.c_str(), OC[k].icast, OC[k].wsize);
    if ( ic != k ) {
      sprintf(MSGERR1, "Invalid cast for column '%s'", OC[k].name.c_str());
      sprintf(MSGERR2, "ICAST=%d", OC[k].icast);
      errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
    }
  }

  // ---- pass 2: build and write one output column at a time
  std::vector<char> out ;
  for(int k=0; k < NCOL_OUT; k++ ) {
    OUTCOL &o = OC[k];
    size_t w  = (size_t)o.wsize ;
    out.assign( NROW_OUT * w + 1, 0 );
    size_t orow = 0 ;

    for(ifile=0; ifile < NFILE; ifile++ ) {
      int i_in = ( o.carried ) ? IMAP[ifile][o.ivar] : -1 ;
      const char *base = NULL ;
      npzrd::COL *pc = NULL;

      if ( i_in >= 0 ) {
	pc = &F[ifile].col[i_in] ;
	if ( !npzrd::inflate_col(FP[ifile], *pc, buf, 0) || buf.size() < pc->npy_off + (size_t)pc->wsize*pc->nrow ) {
	  sprintf(MSGERR1, "Could not read column '%s'", pc->name.c_str());
	  sprintf(MSGERR2, "%.200s", A->INPFILE[ifile]);
	  errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
	}
	base = &buf[0] + pc->npy_off ;
      }

      const std::vector<int> &rows = ROWS[ifile] ;
      for(size_t n=0; n < rows.size(); n++, orow++ ) {
	char *dst = &out[orow*w] ;
	int  m    = ILIST_OUT[orow] ;

	if ( o.src_stored >= 0 && m >= 0 ) 
	  { npz_put_from_stored(dst, o.icast, o.wsize, o.src_stored, m); }
	else if ( !o.carried )        // appended column and undefined row
	  { npz_put_fill(dst, o.icast, o.wsize, VAL_UNDEF_APPENDCOL); }
	else if ( pc != NULL )        // original value from input table
	  { npz_put_from_src(dst, o.icast, o.wsize, *pc, base + (size_t)rows[n]*pc->wsize); }
	else                          // column is not in this file
	  { npz_put_fill(dst, o.icast, o.wsize, -9); }
      }
    }
    npz_table_putcol(h, k, &out[0], NROW_OUT);
  }
  out.clear();  out.shrink_to_fit();

  // ---- header dictionary (JSON in zip comment)
  if ( A->NHDR > 0 ) {
    std::string json = "{";
    for(int j=0; j < A->NHDR; j++ ) {
      if ( j > 0 ) { json += ", "; }
      json += "\"";  npz_json_escape(json, A->HDR_KEY[j]);  json += "\": \"";
      npz_json_escape(json, A->HDR_VAL[j]);  json += "\"";
    }
    json += "}";
    npz_table_set_header(h, json.c_str());
  }

  for(ifile=0; ifile < NFILE; ifile++ ) { fclose(FP[ifile]); }

  printf("   Close NPZ table file: %s \n", A->OUTFILE);  fflush(stdout);
  npz_table_close(h);

  appendcol_print_stats(NROW_INP, (int)NROW_OUT, NMATCH, NUNMATCHED);
  appendcol_hash_free(&H);

  return;

} // end SNTABLE_APPENDCOL_EXEC_NPZ
