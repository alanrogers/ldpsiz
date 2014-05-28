/**
 * @file assign.c
 * @author Alan R. Rogers
 * @brief Functions for objects of class Assignment, which handles
 * assignment statements in input files. These functions do not parse
 * the input files---that job is handled elsewhere, because syntax
 * differs in different files.
 *
 * @copyright Copyright (c) 2014, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "assign.h"
#include "misc.h"
#include "tokenizer.h"

/**
 * A linked list, each item representing the name and value of a
 * variable whose value is assigned in an initialization file. 
 */
struct Assignment {
    char       *key;         /**< name of variable */
    char       *value;       /**< value of variable */
    struct Assignment *next; /**< next Assignment in linked list */
};

/** Insert new assignment. */
Assignment *Assignment_new(Assignment * head, const char *key,
                           const char *value) {
    if(head) {
        int         cmp = strcmp(head->key, key);

        if(cmp == 0) {
            free(head->value);
            head->value = strdup(value);
            return head;
        } else if(cmp < 0) {
            head->next = Assignment_new(head->next, key, value);
            return head;
        }
    }
    /* else head==NULL or cmp > 0 */
    Assignment *new = malloc(sizeof(Assignment));

    checkmem(new, __FILE__, __LINE__);

    new->key = strdup(key);
    new->value = strdup(value);
    new->next = head;

    return new;
}

/** Free list of Assignment objects. */
void Assignment_free(Assignment * head) {
    if(head == NULL)
        return;

    Assignment_free(head->next);
    free(head->key);
    free(head->value);
    free(head);
}

/**
 * Return value corresponding with given key, or NULL if no such key
 * is found.
 */
const char *Assignment_value(const Assignment * head, const char *key) {
    if(head == NULL)
        return NULL;

    int         cmp = strcmp(head->key, key);

    if(cmp == 0)
        return head->value;
    else if(cmp < 0)
        return Assignment_value(head->next, key);
    return NULL;
}

void Assignment_print(const Assignment * head, FILE * fp) {
    if(head == NULL)
        return;

    fprintf(fp, "%-18s = %s\n", head->key, head->value);
    Assignment_print(head->next, fp);
}

/**
 * Initialize *ptr do value corresponding to given key.
 *
 * @param[in] a A structure containing the assigned values.
 *
 * @param[in] mandatory If mandatory!=0, the function aborts if "key"
 * is not found.
 *
 * @param[out] ptr Points to a double variable, which will be
 * initialized. If variable "key" has been assigned, then its value
 * will be written into *ptr. Otherwise, *ptr is unaltered.
 *
 * @returns 1 if variable was initialized, 0 otherwise.
 */
int Assignment_setDbl(const Assignment * a, const char *key, double *ptr,
                      int mandatory) {
    if(!a)
        return 0;
    const char *value = Assignment_value(a, key);

    if(!value) {                /* key not found */
        if(mandatory)
            eprintf("ERR@%s:%d: \"%s\" not assigned",
                    __FILE__, __LINE__, key);
        return 0;
    }
    assert(strlen(value) > 0);
    *ptr = strtod(value, NULL);
    return 1;
}

/**
 * Initialize *ptr do value corresponding to given key.
 *
 * @param[in] a A structure containing the assigned values.
 *
 * @param[in] mandatory If mandatory!=0, the function aborts if "key"
 * is not found.
 *
 * @param[out] ptr Points to a variable of type "long int", which will
 * be initialized. If variable "key" has been assigned, then its value
 * will be written into *ptr. Otherwise, *ptr is unaltered.
 *
 * @returns 1 if variable was initialized, 0 otherwise.
 */
int Assignment_setLong(const Assignment * a, const char *key, long *ptr,
                       int mandatory) {
    if(!a)
        return 0;
    const char *value = Assignment_value(a, key);

    if(!value) {                /* key not found */
        if(mandatory)
            eprintf("ERR@%s:%d: \"%s\" not assigned",
                    __FILE__, __LINE__, key);
        return 0;
    }
    assert(strlen(value) > 0);
    *ptr = strtol(value, NULL, 10);
    return 1;
}

/**
 * Initialize *ptr do value corresponding to given key.
 *
 * @param[in] a A structure containing the assigned values.
 *
 * @param[in] mandatory If mandatory!=0, the function aborts if "key"
 * is not found.
 *
 * @param[out] ptr Points to a variable of type "int", which will be
 * initialized. If variable "key" has been assigned, then its value
 * will be written into *ptr. Otherwise, *ptr is unaltered.
 *
 * @returns 1 if variable was initialized, 0 otherwise.
 */
int Assignment_setInt(const Assignment * a, const char *key, int *ptr,
                      int mandatory) {
    if(!a)
        return 0;
    const char *value = Assignment_value(a, key);

    if(!value) {                /* key not found */
        if(mandatory)
            eprintf("ERR@%s:%d: \"%s\" not assigned",
                    __FILE__, __LINE__, key);
        return 0;
    }
    assert(strlen(value) > 0);
    *ptr = (int) strtol(value, NULL, 10);
    return 1;
}

/**
 * Initialize *ptr do value corresponding to given key.
 *
 * @param[in] a A structure containing the assigned values.
 *
 * @param[in] mandatory If mandatory!=0, the function aborts if "key"
 * is not found.
 *
 * @param[out] ptr Points to a variable of type "unsigned int", which
 * will be initialized. If variable "key" has been assigned, then its
 * value will be written into *ptr. Otherwise, *ptr is unaltered.
 *
 * @returns 1 if variable was initialized, 0 otherwise.
 */
unsigned Assignment_setUnsignedInt(const Assignment * a, const char *key,
                                   unsigned *ptr, int mandatory) {
    if(!a)
        return 0;
    const char *value = Assignment_value(a, key);

    if(!value) {                /* key not found */
        if(mandatory)
            eprintf("ERR@%s:%d: \"%s\" not assigned",
                    __FILE__, __LINE__, key);
        return 0;
    }
    assert(strlen(value) > 0);
    *ptr = (unsigned) strtoul(value, NULL, 10);
    return 1;
}

/**
 * Initialize *ptr do value corresponding to given key.
 *
 * @param[in] a A structure containing the assigned values.
 *
 * @param[in] mandatory If mandatory!=0, the function aborts if "key"
 * is not found.
 *
 * @param[out] ptr Points to a character array containing "size"
 * bytes. If the variable "key" has been assigned, its value will be
 * copied into "ptr". Otherwise, *ptr is unaltered.
 *
 * @returns 1 if variable was initialized, 0 otherwise.
 */
int Assignment_setString(const Assignment * a, const char *key, char *ptr,
                         int size, int mandatory) {
    if(!a)
        return 0;
    const char *value = Assignment_value(a, key);

    if(!value) {                /* key not found */
        if(mandatory)
            eprintf("ERR@%s:%d: \"%s\" not assigned",
                    __FILE__, __LINE__, key);
        return 0;
    }
    assert(strlen(value) > 0);
    snprintf(ptr, size, "%s", value);
    return 1;
}

/**
 * Read assignment statements in file produced by eld.
 *
 * @param[in] ifp Pointer to input file.
 *
 * @returns Pointer to an Assignment object, which contains all the
 * assignments found.
 */
Assignment *Assignment_readEld(FILE * ifp) {
    char        buff[500];
    int         ntokens = -1;
    Tokenizer  *tkz = Tokenizer_new(50);
    Assignment *a = NULL;

    while(1) {
        if(fgets(buff, sizeof(buff), ifp) == NULL)
            break;

        if(!strchr(buff, '\n') && !feof(ifp))
            eprintf("ERR@%s:%d: input buffer overflow."
                    " buff size: %d\n", __FILE__, __LINE__, sizeof(buff));

        if(strstr(buff, "Estimates of sigdsq") != NULL) /* end of header */
            break;

        Tokenizer_split(tkz, buff, "=");    /* tokenize */
        ntokens = Tokenizer_strip(tkz, " \t\n#");   /* strip extraneous */
        if(ntokens != 2)
            continue;           /* not an assignment, so skip this line */

        /*
         * Read assignment statements, store in Assignment structure.
         */
        a = Assignment_new(a,
                           Tokenizer_token(tkz, 0), Tokenizer_token(tkz, 1));
    }

    Tokenizer_free(tkz);
    tkz = NULL;
    return a;
}

/* Local Variables: */
/* mode: c */
/* End: */
