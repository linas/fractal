//
// FILE:
// top.h
//
// FUNCTION:
// Utility macro to keep a top-ten list
//
// HISTORY:
// Linas Vepstas January 1997

#ifndef __LAG_TOP_H__
#define __LAG_TOP_H__

#define UPDATE_TOP_TEN(count,node) {				\
   if (count > top_ten_count[LAG_TOP_TEN-1]) {			\
      for (int i=0; i<LAG_TOP_TEN; i++) {			\
         if (count > top_ten_count[i]) {			\
								\
            /* if slot is empty, use it */			\
            if (!topten[i]) {					\
               top_ten_count [i] = count;			\
               topten[i] = node;				\
               break;						\
            }							\
								\
            /* search for matching id */			\
            int j=0;						\
            int last = LAG_TOP_TEN -1;				\
            for (j=0; j<LAG_TOP_TEN; j++) {			\
               if (topten[j]) {					\
                  if (topten[j]->id == node->id) last = j;	\
               }						\
            }							\
                 						\
            /* this avoids duplicates in the list */		\
            for (j=last; j>i; j--) {				\
               topten[j] = topten[j-1];				\
               top_ten_count[j] = top_ten_count[j-1];		\
            }							\
            top_ten_count[i] = count;				\
            topten[i] = node;					\
            break;						\
         }							\
      }								\
   }								\
}

#endif // __LAG_TOP_H__
