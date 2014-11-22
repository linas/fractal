
#include <fcntl.h>
#include <sys/types.h>
#include <sys/param.h>
#include <sys/filsys.h>
#include <sys/ino.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/dir.h>

/*
 * Attempt to write a UNIX undelete function
 *
 * Linas Vepstas march 1993
 */

main (argc, argv)
int argc;
char *argv[];
{
   char *dev;
   int fd;
   int retval;
   char buff[2048];
   char superblock[2048];
   filsys_t *super;
   int i, j, b, c, l;
   int skip_blocks, inodes_per_block, ino;
   struct dinode *inod;

   dev = "/dev/hd7";

   if (argc !=2) {
      printf ("Usage: %s <inode> \n", argv[0]);
      exit (1);
   }

   ino = atoi (argv[1]);

   fd = open (dev, O_RDONLY);
   if (fd <0) {
      perror ("couldn't open");
      exit (1);
   }

   retval = read (fd, buff, 2048);
   if (retval <0) {
      perror("couldn't read 0 block");
      exit (1);
   }

   retval = read (fd, superblock, 2048);
   if (retval <0) {
      perror("couldn't read superblock");
      exit (1);
   }

   super = (filsys_t *) superblock;

   printf ("Found the superblock \n");
   printf ("file system name %s \n", super->s_fname);
   printf ("volume name %s \n", super->s_fpack);
   printf ("block size (in bytes) %d \n", super->s_bsize);
   printf ("file system size (in blocks) %d \n", super->s_fsize);
   printf ("inode-list size (in blocks) %d \n", super->s_isize);
   printf ("num blocks per cylinder %d \n", super->s_cyl);
   printf ("block interleave (skip) factor %d \n", super->s_skip);
   printf ("num slots in free block list %d \n", super->s_nicfree);
   printf ("num slots in free inode list %d \n", super->s_nicino);
   printf ("num slots in fragment table %d \n", super->s_nicfrag);
   printf ("byte offset to free block list %d \n", super->s_sicfree);
   printf ("byte offset to free inode list %d \n", super->s_sicino);
   printf ("byte offset to fragment table %d \n", super->s_sicfrag);
   printf ("total free blocks %d \n", super->s_tfree);
   printf ("total free inodes %d \n", super->s_tinode);
   printf ("\n");
   printf ("\n");

   printf ("num free block slots %d\n", NICFREE);
   for (i=0; i<NICFREE; i++) {
      printf ("free block %d at 0x%x or 0x%x or %d \n", i, super->s_free[i],
            FREEblk(super, i), FREEblk(super, i));

   }

   printf ("\n");
   printf ("\n");

   printf ("num free inode slots %d\n", NICINOD);
   for (i=0; i<NICINOD; i++) {
      printf ("free inode %d at 0x%x or 0x%x or %d \n", i, super->s_inode[i],
            FREEino(super, i), FREEino(super, i));

   }

   /* compute which  block the desired inode sits in */
   inodes_per_block = BSIZE / 64;
   skip_blocks = ino / inodes_per_block;

#ifdef MANUAL
   for (b=0; b<skip_blocks; b++) {
      retval = read (fd, buff, 2048);
      if (retval <0) {
         fprintf (stderr, " block %d ", b);
         perror("couldn't read ");
         exit (1);
      }
      printf ("successfully read block %d +2 \n", b);
   }
#endif 
   
   b = skip_blocks +2;
   retval = lseek (fd, b * BSIZE, SEEK_SET);
   if (retval < 0) {
      perror ("seek for inode block failed");
      exit (1);
   }

   retval = read (fd, buff, 2048);
   if (retval <0) {
      fprintf (stderr, " block %d ", b);
      perror("couldn't read ");
      exit (1);
   }
   printf ("\n");
   printf ("\n");
   printf ("successfully read block %d containing desired inode \n", b);

   printf ("searching for inode %d ", ino);
   ino -= skip_blocks * inodes_per_block;
   ino --;
   printf ("the %d th node in this block \n", ino);

   inod = (struct dinode *) buff;
   inod = &inod[ino];

   printf (" file mode is 0x%x \n", inod -> di_mode);
   printf (" num links to file %d \n", inod -> di_nlink);
   printf (" owner uid %d \n", inod -> di_uid);
   printf (" owner gid %d \n", inod -> di_gid);
   printf (" size in bytes %d \n", inod -> di_size);
   printf (" last accessed on %s \n", ctime (&(inod -> di_atime)));
   printf (" last modified on %s \n", ctime (&(inod -> di_mtime)));
   printf (" created on %s \n", ctime (&(inod -> di_ctime)));

   printf ("\n");
   printf ("disk block addresses are: \n");
   for (i=0; i<39; i += 3) {
      j = inod->di_addr [i];
      j <<= 8;
      j |= inod->di_addr [i+1];
      j <<= 8;
      j |= inod->di_addr [i+2];
      printf ("0x%x or %d \n", j, j);
   }
      

   i=0;
   j = inod->di_addr [i];
   j <<= 8;
   j |= inod->di_addr [i+1];
   j <<= 8;
   j |= inod->di_addr [i+2];

   for (i=0; i<NICFREE; i++) {
      printf ("\n");
      printf ("\n");
      b = FREEblk(super, i);
      printf ("free block %d at 0x%x or %d \n", i, b, b);

      retval = lseek (fd, b * BSIZE, SEEK_SET);
      if (retval < 0) {
         perror ("seek for file block failed");
         exit (1);
      }

      retval = read (fd, buff, 2048);
      if (retval <0) {
         fprintf (stderr, " block %d ", b);
         perror("couldn't read ");
         exit (1);
      }

      printf ("successfully read block %d \n", b);
      printf ("dumping block %d \n", b);
      printf ("\n");


      for (l=0; l<32; l++) {
         printf ("%d  ",64*l );
         for (c=0; c<64; c++) {
            printf ("%c", buff[64*l+c]);
         }
         printf ("\n");

      }

      /*
      for (l=0; l<128; l++) {
         struct direct *d;
   
         d = (struct direct *) &buff[l*16];
         printf ("entry %d has inode %d name %s \n", l, d->d_ino, d->d_name);
      }
      */
   }
}


