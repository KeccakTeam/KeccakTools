/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is based on genKAT.c by NIST.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#include "Keccak.h"

#define MAX_MARKER_LEN      50
#define SUBMITTER_INFO_LEN  128

typedef enum { KAT_SUCCESS = 0, KAT_FILE_OPEN_ERROR = 1, KAT_HEADER_ERROR = 2, KAT_DATA_ERROR = 3, KAT_HASH_ERROR = 4 } STATUS_CODES;
typedef unsigned char BitSequence;

#define SqueezingOutputLength 4096

STATUS_CODES genShortMsg(int hashbitlen);
int     FindMarker(FILE *infile, const char *marker);
int     ReadHex(FILE *infile, BitSequence *A, int Length, char *str);
void    fprintBstr(FILE *fp, char *S, BitSequence *A, int L);


void genKATShortMsg_main()
{
    genShortMsg(0);
    genShortMsg(224);
    genShortMsg(256);
    genShortMsg(384);
    genShortMsg(512);
}

void fromNISTConventionToInternalConvention(BitSequence *Msg, int msglen)
{
    if ((msglen % 8) != 0) {
        Msg[msglen/8] >>= 8 - (msglen % 8);
    }
}

STATUS_CODES genShortMsg(int hashbitlen)
{
    char        fn[32], line[SUBMITTER_INFO_LEN];
    int         msglen, msgbytelen, done;
    BitSequence Msg[256], MD[64];
    BitSequence Squeezed[SqueezingOutputLength/8];
    FILE        *fp_in, *fp_out;
    
    if ( (fp_in = fopen("ShortMsgKAT.txt", "r")) == NULL ) {
        printf("Couldn't open <ShortMsgKAT.txt> for read\n");
        return KAT_FILE_OPEN_ERROR;
    }
    
    sprintf(fn, "ShortMsgKAT_%d.txt", hashbitlen);
    if ( (fp_out = fopen(fn, "w")) == NULL ) {
        printf("Couldn't open <%s> for write\n", fn);
        return KAT_FILE_OPEN_ERROR;
    }
    fprintf(fp_out, "# %s\n", fn);
    if ( FindMarker(fp_in, "# Algorithm Name:") ) {
        fscanf(fp_in, "%[^\n]\n", line);
        fprintf(fp_out, "# Algorithm Name:%s\n", line);
    }
    else {
        printf("genShortMsg: Couldn't read Algorithm Name\n");
        return KAT_HEADER_ERROR;
    }
    if ( FindMarker(fp_in, "# Principal Submitter:") ) {
        fscanf(fp_in, "%[^\n]\n", line);
        fprintf(fp_out, "# Principal Submitter:%s\n", line);
    }
    else {
        printf("genShortMsg: Couldn't read Principal Submitter\n");
        return KAT_HEADER_ERROR;
    }
    
    done = 0;
    do {
        if ( FindMarker(fp_in, "Len = ") )
            fscanf(fp_in, "%d", &msglen);
        else {
            done = 1;
            break;
        }
        msgbytelen = (msglen+7)/8;

        if ( !ReadHex(fp_in, Msg, msgbytelen, "Msg = ") ) {
            printf("ERROR: unable to read 'Msg' from <ShortMsgKAT.txt>\n");
            return KAT_DATA_ERROR;
        }
        fprintf(fp_out, "\nLen = %d\n", msglen);
        fprintBstr(fp_out, "Msg = ", Msg, msgbytelen);
        Keccak keccak(  (hashbitlen <= 256) ? 1024 :  512, 
                        (hashbitlen <= 256) ?  576 : 1088,
                        hashbitlen/8);
        fromNISTConventionToInternalConvention(Msg, msglen);
        keccak.absorb(Msg, msglen);
        keccak.pad();
        if (hashbitlen > 0)
            keccak.squeeze(MD, hashbitlen);
        else
            keccak.squeeze(Squeezed, SqueezingOutputLength);
        if (hashbitlen > 0)
            fprintBstr(fp_out, "MD = ", MD, hashbitlen/8);
        else
            fprintBstr(fp_out, "Squeezed = ", Squeezed, SqueezingOutputLength/8);
    } while ( !done );
    printf("finished ShortMsgKAT for <%d>\n", hashbitlen);
    
    fclose(fp_in);
    fclose(fp_out);
    
    return KAT_SUCCESS;
}

//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
int
FindMarker(FILE *infile, const char *marker)
{
    char    line[MAX_MARKER_LEN];
    int     i, len;

    len = (int)strlen(marker);
    if ( len > MAX_MARKER_LEN-1 )
        len = MAX_MARKER_LEN-1;

    for ( i=0; i<len; i++ )
        if ( (line[i] = fgetc(infile)) == EOF )
            return 0;
    line[len] = '\0';

    while ( 1 ) {
        if ( !strncmp(line, marker, len) )
            return 1;

        for ( i=0; i<len-1; i++ )
            line[i] = line[i+1];
        if ( (line[len-1] = fgetc(infile)) == EOF )
            return 0;
        line[len] = '\0';
    }

    // shouldn't get here
    return 0;
}

//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
int
ReadHex(FILE *infile, BitSequence *A, int Length, char *str)
{
    int         i, ch, started;
    BitSequence ich;

    if ( Length == 0 ) {
        A[0] = 0x00;
        return 1;
    }
    memset(A, 0x00, Length);
    started = 0;
    if ( FindMarker(infile, str) )
        while ( (ch = fgetc(infile)) != EOF ) {
            if ( !isxdigit(ch) ) {
                if ( !started ) {
                    if ( ch == '\n' )
                        break;
                    else
                        continue;
                }
                else
                    break;
            }
            started = 1;
            if ( (ch >= '0') && (ch <= '9') )
                ich = ch - '0';
            else if ( (ch >= 'A') && (ch <= 'F') )
                ich = ch - 'A' + 10;
            else if ( (ch >= 'a') && (ch <= 'f') )
                ich = ch - 'a' + 10;

            for ( i=0; i<Length-1; i++ )
                A[i] = (A[i] << 4) | (A[i+1] >> 4);
            A[Length-1] = (A[Length-1] << 4) | ich;
        }
    else
        return 0;

    return 1;
}

void
fprintBstr(FILE *fp, char *S, BitSequence *A, int L)
{
    int     i;

    fprintf(fp, "%s", S);

    for ( i=0; i<L; i++ )
        fprintf(fp, "%02X", A[i]);

    if ( L == 0 )
        fprintf(fp, "00");

    fprintf(fp, "\n");
}
