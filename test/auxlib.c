// SPDX-License-Identifier: BSD-3-Clause
/**
 * A tiny lib of string functions that C misses
 * Created on Sun Dec 29 21:34:24 2024
 * 
 * @author: hilton
 * 
 **/
 
#include <stdlib.h> /* calloc(), realloc(), free() */
#include <string.h> /* strcpy(), strlen() */
#include <regex.h>  /* regex_t, regmatch_t, regcomp(), regexec() , 
					   regfree() */

#include "auxlib.h"

/**
 * Strips the space characters of a string to the left and right
 * 
 * @param[in] input The string to be deblanked
 * 
 * \retval A new string with left and right spaces eliminated 
 * In case of sucess.
 * 
 * \retval The null char pointer 
 * In case of failure.
 **/
  
char* 
deblank(const char* input) 
{
	int start = 0;
	int stop = 0;
	size_t length = 0;
	size_t full_length = 0;
	char* result = (char*)NULL;
	
	if (!input) return NULL;
	
	full_length = strlen(input);
	if (0 == full_length) {
		return (char*)calloc(1, sizeof(char*));
	}
	
	for (start = 0; *(input + start) == ' '; start++)
		;
	
	for (stop = full_length - 1; stop > 0; stop--) {
		if (*(input + stop) != ' ') break;
	} 
	
	length = full_length - start + 1; 
	length -= full_length - stop;
	
	result = (char*)calloc(length + 1, sizeof(char));
	strncpy(result, input + start, length);
	
	/* Normal function termination */
	return result;
}


/**
 * Checks if a string contains only a number.
 * Spaces to the left and right are ignored.
 * 
 * @param[in] input The string to be verified.
 * 
 * \retval 1 
 * In case of sucess.
 * 
 * \retval 0
 * If the string doesn't contain only a number
 **/

int 
is_numeric(const char* input) 
{
	int rv;
	int result;
	char* deb = (char*)NULL;
	regex_t regex;
	regmatch_t pmatch[3];         // Up to 3 sub-expressions
	
	rv = regcomp(&regex, NUMERIC_MASK, REG_EXTENDED | REG_ICASE | REG_NEWLINE);
	if (rv) {
		return 0;
	}

	/* printf("input '%s'\n", input); */
	deb = deblank(input);
	
	/* printf("deb '%s'\n", deb); */
	rv = regexec(&regex, deb, ARRAY_SIZE (pmatch), pmatch, 0);
	
	if (rv) {
		return 0;
	}
	
	result = (int)
		(strlen(deb) == (size_t)(pmatch[0].rm_eo - pmatch[0].rm_so));
	
	regfree(&regex);
	free(deb);
	
	/* Normal function termination */
	return result; 
}

/**
 * Splits a string in several others delimited by a separator.
 * 
 * @param[in] text The string containing the text to be splitted.
 * @param[in] sep  A delimiter in the input string.
 * @param[out] fields A char array pointing tho the strings splitted.
 * 
 * \retval The number of substrings found. 
 * In case of sucess. If the string doesn't contain a delimiter, it is 
 * copied to the output, and 1 is returned.
 * 
 * \retval 0
 * In case of some failure. 
 * 
 * \note Uses `strtoken()` and is therefore not reentrant.
 **/

size_t 
split_fields(const char* text, const char* sep, char*** fields) 
{
	char*  buf    = (char*)NULL;
	char*  token  = (char*)NULL;
	char**  work  = (char**)NULL;
	size_t result = 0;
	size_t n_field = 0;
	
	buf = (char*)calloc(strlen(text) + 1, sizeof(char));
	strcpy(buf, text);
	
	token = strtok(buf, sep);
    /* TODO if separator not found, leave */
    if (!token) {
		/* TODO report reading problem */
        return 0;
	}
	// printf("1st token '%s'\n", token);
	
	work = (char**)calloc(1024, sizeof(char*));

    /* TODO get all elements */
	while (token) {
		work[n_field] = 
			(char*)calloc(strlen(token) + 1, sizeof(char*));
		strcpy(work[n_field], token);
		
		n_field++;
		token = strtok(NULL, sep);
		// printf("token %ld '%s'\n", n_field, token);
	}

	result = n_field;
	
	free(buf);
	work = (char**)realloc(work, result * sizeof(char*));
	
	*fields = work;
	
	/* Normal function termination */
	return result;
}

/**
 * Free an array of strings, just as the ones returned by 
 * `split_fields()`
 * 
 * @param[in] n_fields The number of elements in the string array.
 * @param[out] fields A string array to be freed.
 * 
 * \retval 0. 
 * In case of sucess. 
 * 
 * \retval 1
 * In case of some failure, like zero number of elements, or 
 * if the char array is the null pointer
 **/

int 
free_fields(size_t n_fields, char*** fields) 
{
	size_t ind;
	char** work;
	
	if ((n_fields == 0) || (! fields) || (! *fields)) {
		return 1;
	}
	
	work = *fields;
	
	for (ind = 0; ind < n_fields; ind++) {
		if (! work[ind]) {
			continue;
		}
		free(work[ind]);
	}
	
	free(work);
	
	return 0;
}

/**
 * Search a text in an array of strings, just as the ones returned by 
 * `split_fields()`
 * 
 * @param[in] texts The array of strings.
 * @param[in] n_texts The number of texts in the string array.
 * @param[int] key The text to be freed. Should not be the empty string.
 * 
 * \retval The zero-based position where the key is. 
 * In case of sucess. 
 * 
 * \retval -1
 * In case text is not found or some failure
 * of some failure, like zero number of elements, or 
 * if the char array is the null pointer
 **/

int 
find_text(char* const* texts, size_t n_texts, const char* key) 
{
	size_t ind;
	
	/* if there's not where to find, return nothing was found */
	if ((NULL == texts) || (0 == n_texts)) {
		return -1;
	}
	
	/* if there's not what to find, return nothing was found */
	if ( (NULL == key) || (0 == strlen(key)) ) {
		return -1;
	}
	
	for (ind = 0; ind < n_texts; ind++) {
		if (!strcmp(texts[ind], key)) {
			return ind;
		}
	}
	
	return -1;
}
