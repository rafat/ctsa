// SPDX-License-Identifier: BSD-3-Clause
/*
 * Tiny library to read CSV files in C
 * Created on Fri Dec 27 19:37:13 2024
 * 
 * @author: hilton
 * 
 * OBS C99 at least is recommended 
 */

#include <stdio.h>   /* fopen(), fclose() */
#include <stdlib.h>  /* calloc(), free() */
#include <string.h>  /* strtok(), strlen(), strcnmp() */
#include <stdbool.h> /* bool, true, false */

#include "csv_in.h" /* NAN, read_csv_column */
#include "auxlib.h" /* deblank(), find_text(), is_numeric(), 
					   split_fields(), free_fields() */

#define QUANTUM 2048

#if FUTURE_VERSION
struct _CSV {
    
};

/* TODO struct object-like CSV with open that opens the file and reads the 1st line */
/* TODO holds the field names, returns some columns, closes the file */

#endif /* 0 */

/**
 * Create a double vector from a CSV file name, given its column.
 * 
 * @param[in]  csv_name Name of input file.
 * @param[in]  col_name Name of column to extract.
 * @param[out] column `double` vector that contains the column data.
 * 
 * \retval The number of rows in the output vector.
 * In case of success. 
 * 
 * \retval 0
 * In case of failure. 
 **/
 
size_t 
read_csv_column(const char* csv_name, const char* col_name, double** column) 
{
	int rv;
    char buf[1024] = "";
    char* pbuf = (char*)NULL;
    char** row     = (char**)NULL;
    char** columns = (char**)NULL;
    double* work = NULL;
    size_t row_no = 0;
    size_t n_fields = 0;
    size_t n_columns = 0;
    size_t col_no = (size_t)-1;
    size_t total_rows = 0;
    size_t result = 0;
    FILE* csv_f = NULL;
    
    /* HINT opening the input file*/
    csv_f = fopen(csv_name, "r");
    if (NULL == csv_f) {
		/* TODO report opening problem */
		printf("ERROR Could not open file '%s' for input\n", csv_name);
		
        return (size_t)0;
    }
    
    /* HINT getting the header */
    pbuf = fgets(buf, sizeof(buf) - 1, csv_f);
    if (NULL == pbuf) {
		/* TODO report reading problem */
		printf("ERROR Could not read 1st line\n");
		
        return (size_t)0;
    }
	
	n_columns = split_fields(buf, ",", &columns);
	
	rv = find_text(columns, n_columns, col_name);
	
    /* TODO if the suggested column name is not in the header, leave */
    if (-1 == rv) {
		/* TODO report reading problem */
		printf("ERROR Column '%s' was not found in header\n", col_name);
		
		puts("Columns found\n");
		for (col_no = 0; col_no < n_columns; col_no++) {
			printf("'%s' ", columns[col_no]);
		}
		
        return (size_t)0;
	}
	
	free_fields(n_columns, &columns);

    
    col_no = (size_t)rv;
    
    /* TODO allocate a number of elements in the workorary area */
    work = (double*)calloc(QUANTUM, sizeof(double));
    total_rows = QUANTUM;
    
    /* TODO for the rest of the file do */
    for (pbuf = fgets(buf, sizeof(buf) - 1, csv_f); 
		NULL != pbuf;
		pbuf = fgets(buf, sizeof(buf) - 1, csv_f)) { 
			 
		/* TODO break the CSV row in comma-separated values */
		n_fields = split_fields(buf, ",", &row);
		
		/* TODO if work size is too small, realloc the work array with QUANTUM more */
		if (total_rows <= (row_no + 1)) {
			work = realloc(work, total_rows + QUANTUM);
			total_rows += QUANTUM;
		}
		
		/* TODO if the column number > number of fields, set NAN to it */
		if (col_no > (n_fields - 1)) {
			work[row_no] = NAN;
		}
		
		/* TODO if there's no element in the suggested position, set NAN to it */
		if ( ! strlen(row[col_no]) ) {
			work[row_no] = NAN;
		}
		
		/* TODO if there's no valid number in the suggested position, set NAN to it */
		if ( ! is_numeric(row[col_no]) ) {
			work[row_no] = NAN;
		}

		/* TODO add the number to the work area*/
		work[row_no] = strtod(row[col_no], NULL);

		row_no++;
		
		free_fields(n_fields, &row);
	}
    
    /* Reallocate the result area, using the number of elements */
    work = (double*)realloc(work, row_no*sizeof(double));
    
    /* Copy the work area to the result area */
    *column = work;
    
    result = row_no;
    
    /* Close the input file */
    fclose(csv_f);
    
    /* Normal function termination */
    return result;
}
