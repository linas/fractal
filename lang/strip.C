//
// FILE:
// strip.C
//
// FUNCTION:
// strip mail headers
//
// HISTORY:
// January 1997 Linas Vepstas

#include <ctype.h>
#include <stdio.h>
#include <string.h>


main () {

   FILE * fh = stdin;
   
   char buff[5000];

   while (!feof(fh)) {
      fgets (buff, 5000, fh);
      buff[4999] = 0x0;

      // reject the mail headers
      if (!strncmp (buff, "From", 4)) continue;
      if (strstr (buff, "From:")) continue;
      if (strstr (buff, "Date:")) continue;
      if (strstr (buff, "Subject:")) continue;
      if (strstr (buff, "SUBJECT:")) continue;
      if (strstr (buff, "To:")) continue;
      if (strstr (buff, "Re:")) continue;
      if (strstr (buff, "Status:")) continue;
      if (strstr (buff, "MIME-Version:")) continue;
      if (strstr (buff, "MIME-version:")) continue;
      if (strstr (buff, "Content-Type:")) continue;
      if (strstr (buff, "Return-Path:")) continue;
      if (strstr (buff, "Received:")) continue;
      if (strstr (buff, "In-Reply-To:")) continue;
      if (strstr (buff, "Message-ID:")) continue;
      if (strstr (buff, "http:")) continue;
      if (strstr (buff, "www")) continue;
      if (strstr (buff, "Newsgroups:")) continue;
      if (strstr (buff, "Content-Type:")) continue;
      if (strstr (buff, "Content-Transfer-Encoding:")) continue;
      if (strstr (buff, "Content-transfer-encoding:")) continue;
      if (strstr (buff, "Appearently-To:")) continue;

      if (strstr (buff, "Message-Id:")) continue;
      if (strstr (buff, "X-Sender:")) continue;
      if (strstr (buff, "X-Mailer:")) continue;
      if (strstr (buff, "Mime-Version:")) continue;
      if (strstr (buff, "CC:")) continue;
      if (strstr (buff, "Cc:")) continue;
      if (strstr (buff, "cc:")) continue;
      if (strstr (buff, "X-Sun-Charset:")) continue;
      if (strstr (buff, "Sender:")) continue;
      if (strstr (buff, "Orgnaization:")) continue;
      if (strstr (buff, "References:")) continue;
      if (strstr (buff, "X-Loop:")) continue;
      if (strstr (buff, "In-reply-to")) continue;
      if (strstr (buff, "X-UIDL:")) continue;
      if (strstr (buff, "X-MSMail-Priority")) continue;
      if (strstr (buff, "X-Priority")) continue;

      if (strstr (buff, "teleportal")) continue;
/*
      if (strstr (buff, "")) continue;
      if (strstr (buff, "")) continue;
      if (strstr (buff, "")) continue;
      if (strstr (buff, "")) continue;
*/
      // if (strstr (buff, "")) continue;

      printf ("%s", buff);
   }

}

// ========================== END OF FILE ==================
