// Created Apr 2025
// Functions to read/write npz file used by python.

#define MXTABLE_NPZ  10
#define ISFILE_NPZ_COV        1
#define ISFILE_NPZ_FACTORIZE  2

int NTABLE_NPZ ;

#ifdef __cplusplus
extern"C" {
#endif

  int read_npz_covmat(char *npz_file, double *array1d);
  int read_npz_factorized(char *npz_file, double **diag, double **U, int *K);
  int get_npz_filetype(char *npz_file);

  // ---- table writer (buffered, spilled to compressed temp files, zipped at close)
  int  npz_table_create (char *fileName);
  int  npz_table_addcol (int h, char *varName, int icast, int strlen_max);
  void npz_table_append (int h, int icol, void *ptr);
  void npz_table_endrow (int h);     // call once per row; triggers flush every NFLUSH rows
  void npz_table_close  (int h);

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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

//#include "sntools_npz.h"
#include "cnpy.h"


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
    NPZTab() : nrow(0), nbuf(0), nflush(100), level(6), used(false) {}
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

  std::vector<char> npy_header_for(int icast, size_t nrow, int wsize) {
    std::vector<size_t> s1 = { nrow };
    std::vector<size_t> s2 = { nrow, (size_t)wsize };
    switch (icast) {
      case NPZ_ICAST_D: return cnpy::create_npy_header<double>(s1);
      case NPZ_ICAST_F: return cnpy::create_npy_header<float>(s1);
      case NPZ_ICAST_I: return cnpy::create_npy_header<int>(s1);
      case NPZ_ICAST_S: return cnpy::create_npy_header<short>(s1);
      case NPZ_ICAST_L: return cnpy::create_npy_header<long long>(s1);
      default:          return cnpy::create_npy_header<char>(s2);   // fixed-width strings
    }
  }

  void put64(std::vector<char> &v, uint64_t x) { for(int i=0;i<8;i++) v.push_back((char)((x>>(8*i))&0xFF)); }
} // end anonymous namespace


int npz_table_create(char *fileName) {
  for (int h=0; h < MXTABLE_NPZ; h++) {
    if ( !NPZ_TAB[h].used ) {
      NPZ_TAB[h] = NPZTab();  NPZ_TAB[h].used = true;  NPZ_TAB[h].file = fileName;
      const char *e;
      if ( (e = getenv("SNANA_NPZ_NFLUSH")) && atoi(e) > 0 ) { NPZ_TAB[h].nflush = atoi(e); }
      if ( (e = getenv("SNANA_NPZ_LEVEL"))  && atoi(e) >= 0 && atoi(e) <= 9 ) { NPZ_TAB[h].level = atoi(e); }
      return h;
    }
  }
  return -1;
}

int npz_table_addcol(int h, char *varName, int icast, int strlen_max) {
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
  NPZCol &c = NPZ_TAB[h].col[icol];
  size_t n0 = c.buf.size();
  c.buf.resize(n0 + c.wsize, 0);
  if ( c.icast == NPZ_ICAST_C ) { memcpy(&c.buf[n0], ptr, strnlen((char*)ptr, c.wsize)); }
  else                          { memcpy(&c.buf[n0], ptr, c.wsize); }
}

void npz_table_endrow(int h) {
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
    put16(tail, 0);

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
}
