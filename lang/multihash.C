//
// FILE:
// multihash.C
//
// FUNCTION:
// Define multiple word tuplets
//
// HISTORY:
// Linas Vepstas January 1997

#include "multihash.h"

// ------------------------
#define lagGenericWordTable lagWordPairTable
#define LAG_TWO_WORD
#define LAG_WORD_TUPLE 2
#include "pairhash.C"
#undef lagGenericWordTable 
#undef LAG_TWO_WORD
#undef LAG_WORD_TUPLE

// ------------------------
#define lagGenericWordTable lagWordTripleTable
#define LAG_THREE_WORD
#define LAG_WORD_TUPLE 3
#include "pairhash.C"
#undef lagGenericWordTable 
#undef LAG_THREE_WORD
#undef LAG_WORD_TUPLE

// ------------------------
#define lagGenericWordTable lagWordQuadTable
#define LAG_FOUR_WORD
#define LAG_WORD_TUPLE 4
#include "pairhash.C"
#undef lagGenericWordTable 
#undef LAG_FOUR_WORD
#undef LAG_WORD_TUPLE

// ------------------------
#define lagGenericWordTable lagWordQuintTable
#define LAG_FIVE_WORD
#define LAG_WORD_TUPLE 5
#include "pairhash.C"
#undef lagGenericWordTable 
#undef LAG_FIVE_WORD
#undef LAG_WORD_TUPLE

// ------------------------
#define lagGenericWordTable lagWordHexTable
#define LAG_SIX_WORD
#define LAG_WORD_TUPLE 6
#include "pairhash.C"
#undef lagGenericWordTable 
#undef LAG_SIX_WORD
#undef LAG_WORD_TUPLE

// ========================================================
// ------------------------
#define lagGenericWordTable lagWordPairTable
#define lagGenericConcordTable lagConcordPairTable
#include "concord.C"
#undef lagGenericConcordTable 
#undef lagGenericWordTable 

// ------------------------
#define lagGenericWordTable lagWordTripleTable
#define lagGenericConcordTable lagConcordTripleTable
#include "concord.C"
#undef lagGenericConcordTable 
#undef lagGenericWordTable 

// ------------------------
#define lagGenericWordTable lagWordQuadTable
#define lagGenericConcordTable lagConcordQuadTable
#include "concord.C"
#undef lagGenericConcordTable 
#undef lagGenericWordTable 

// ------------------------
#define lagGenericWordTable lagWordQuintTable
#define lagGenericConcordTable lagConcordQuintTable
#include "concord.C"
#undef lagGenericConcordTable 
#undef lagGenericWordTable 

// ------------------------
#define lagGenericWordTable lagWordHexTable
#define lagGenericConcordTable lagConcordHexTable
#include "concord.C"
#undef lagGenericConcordTable 
#undef lagGenericWordTable 

// ================== END OF FILE ==================
