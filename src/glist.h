#ifndef _GLIST_H
#define _GLIST_H

#define LIST(TYPE) struct _list_struct_ ## TYPE * 


#endif


#define LIST_HEADER(LIST_TYPE)\
struct _list_struct_ ## LIST_TYPE {             \
    struct _list_struct_ ## LIST_TYPE * next;   \
    LIST_TYPE val;                              \
};                                              \
\
LIST(LIST_TYPE) _list_add_element_ ## LIST_TYPE (const LIST(LIST_TYPE) oldlist, const LIST_TYPE newelt);\

LIST_HEADER(float)

