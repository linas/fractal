
// FILE:
// multihash.h
//
// FUNCTION:
// Define multiple word tuplets
//
// HISTORY:
// Linas Vepstas January 1997

#ifndef __LAG_MULTI_TABLE_H__
#define __LAG_MULTI_TABLE_H__

// ------------------------
#ifdef __LAG_PAIR_TABLE_H__
#undef __LAG_PAIR_TABLE_H__
#endif // __LAG_PAIR_TABLE_H__

#define lagGenericWordTable lagWordPairTable
#define LAG_TWO_WORD
#define LAG_WORD_TUPLE 2
#include "pairhash.h"
#undef lagGenericWordTable 
#undef LAG_TWO_WORD
#undef LAG_WORD_TUPLE

// ------------------------
#ifdef __LAG_PAIR_TABLE_H__
#undef __LAG_PAIR_TABLE_H__
#endif // __LAG_PAIR_TABLE_H__

#define lagGenericWordTable lagWordTripleTable
#define LAG_THREE_WORD
#define LAG_WORD_TUPLE 3
#include "pairhash.h"
#undef lagGenericWordTable 
#undef LAG_THREE_WORD
#undef LAG_WORD_TUPLE

// ------------------------
#ifdef __LAG_PAIR_TABLE_H__
#undef __LAG_PAIR_TABLE_H__
#endif // __LAG_PAIR_TABLE_H__

#define lagGenericWordTable lagWordQuadTable
#define LAG_FOUR_WORD
#define LAG_WORD_TUPLE 4
#include "pairhash.h"
#undef lagGenericWordTable 
#undef LAG_FOUR_WORD
#undef LAG_WORD_TUPLE


// ------------------------
#ifdef __LAG_PAIR_TABLE_H__
#undef __LAG_PAIR_TABLE_H__
#endif // __LAG_PAIR_TABLE_H__

#define lagGenericWordTable lagWordQuintTable
#define LAG_FIVE_WORD
#define LAG_WORD_TUPLE 5
#include "pairhash.h"
#undef lagGenericWordTable 
#undef LAG_FIVE_WORD
#undef LAG_WORD_TUPLE


// ------------------------
#ifdef __LAG_PAIR_TABLE_H__
#undef __LAG_PAIR_TABLE_H__
#endif // __LAG_PAIR_TABLE_H__

#define lagGenericWordTable lagWordHexTable
#define LAG_SIX_WORD
#define LAG_WORD_TUPLE 6
#include "pairhash.h"
#undef lagGenericWordTable 
#undef LAG_SIX_WORD
#undef LAG_WORD_TUPLE


#endif // __LAG_MULTI_TABLE_H__
