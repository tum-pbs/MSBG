#define BMP_PDS_GEN			0
#define BMP_PDS_MROCTX_RAW		1
//#define BMP_PDS_LROC_CDR	2
#define BMP_PDS_SELENE_DTM		3

typedef struct 
{
  int lines,
      samples,
      sample_bits,
      sample_type,
      sample_dummy,
      sample_dummy_valid;
  unsigned long data_offs;
}
PDS_Info;

int BmpReadPDS(char *fname, 
    	       int pds_typ,
	       float sf,		// scaling factor
               int tmap_typ,		// tone mapping type
    	       BmpBitmap **pBmp );

