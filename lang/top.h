
//
// top.h
//
// Some utility defines
//
// Linas Vepstas January 1997

#define UPDATE_TOP_TEN(count,node) {				\
   if (count > top_ten_count[LAG_TOP_TEN-1]) {			\
      for (i=0; i<LAG_TOP_TEN; i++) {				\
         if (count > top_ten_count[i]) {			\
            for (j=LAG_TOP_TEN-1; j>i; j--) {			\
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
