#ifndef TSDATA_H_
#define TSDATA_H_


#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <ctype.h>

#define BUFFER 16384

#ifdef __cplusplus
extern "C" {
#endif

typedef struct tsdata_set* tsdata_object;

tsdata_object tsdata_read_file(const char *filepath, const char *delimiter, int index, int header);

struct tsdata_set {
	int N;
	double *data;
	double *exog;
	float params[1];
};


#ifdef __cplusplus
}
#endif

#endif /* TSDATA_H_ */