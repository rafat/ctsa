// SPDX-License-Identifier: BSD-3-Clause
/**
 * Header for a tiny lib of string functions that C misses
 * Created on Sun Dec 29 21:34:24 2024
 * 
 * @author: hilton
 * 
 **/

#if !defined(__AUXLIB_H)
#define __AUXLIB_H

#define NUMERIC_MASK "[+-]?[0-9]\\.?[0-9]*(e[+-]?[0-9]+)?"

#define ARRAY_SIZE(x)  ( sizeof(x) / sizeof(x[0]) )

 
char*  deblank(const char* input); 
int    is_numeric(const char* input);
size_t split_fields(const char* text, const char* sep, char*** fields);
int    free_fields(size_t n_fields, char*** fields);
int    find_text(char* const* texts, size_t n_texts, const char* key); 

#endif /* !defined(__AUXLIB_H) */
