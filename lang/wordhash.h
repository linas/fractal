//
// wordhash.C
//
// The lagWordtable class associates a unique ID number with a string.
// The GetWordID() method returns the ID number of the string.
//
// January 1997 Linas Vepstas

#ifndef __LAG_WORD_TABLE_H__
#define __LAG_WORD_TABLE_H__

#define LAG_WORD_TABLE_SIZE 65536

class lagWordTable {
   public:
      lagWordTable (void);
      ~lagWordTable ();
      int GetWordID (char *);
      void Dump (void);

   protected:
      int Hash (char *);

   private:
      class Helper {
         public:
            Helper * next;
            char * theword;
            int id;
            int cnt;
      };
      Helper * table [LAG_WORD_TABLE_SIZE];

      int unused_id;
      int num_entries;
      int num_processed;
      int num_collisions;
};
   
#endif // __LAG_WORD_TABLE_H__

// ================= END OF FILE ==================

