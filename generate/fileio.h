/*
 * fileio.h
 *
 * HISTORY:
 * quick and dirty hack Linas Vepstas October 1989
 * more hacks ever since -- Linas
 */

/*-------------------------------------------------------------------*/
// Read a *.flo greyscale floating-point pixmap.
// The file format is width, height in ascii, then newline
// then floating point data in machine-native byte order.
// One float per pixel, a total of width*height floats.
float* read_flo_file(const char *fname,
                     unsigned int* data_width, unsigned int* data_height);

/*-------------------------------------------------------------------*/
// Read a *.pfm greyscale floating-point pixmap.
// This reads the PFM file format as described in various places.
float* read_pfm_file(const char *fname,
                     unsigned int* data_width, unsigned int* data_height);

/*-------------------------------------------------------------------*/

// Auto-handle above types
float* read_floats(const char *fname, const char** suff,
                   unsigned int* data_width, unsigned int* data_height);

/*-------------------------------------------------------------------*/
// Write a *.flo greyscale floating-point pixmap.
// The file format is width, height in ascii, then newline
// then floating point data in machine-native byte order.
// One float per pixel, a total of width*height floats.
void write_flo_file(const char *fname, float* data,
                    unsigned int data_width, unsigned int data_height);

/*-------------------------------------------------------------------*/
// Write a *.pfm greyscale floating-point pixmap.
// The file format is width, height in ascii, then newline
// then floating point data in machine-native byte order.
// One float per pixel, a total of width*height floats.
void write_pfm_file(const char *fname, float* data,
                    unsigned int data_width, unsigned int data_height);

/*-------------------------------------------------------------------*/
// Auto-handle above types
void write_floats(const char *fname, const char* suff, float* data,
                  unsigned int data_width, unsigned int data_height);

/* --------------- END OF FILE ------------------------- */
