// Created Apr 2025
// Functions to read/write npz file used by python.



#ifdef __cplusplus
extern"C" {
#endif

  int read_npz_covmat(char *npz_file, double *array1d);

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
