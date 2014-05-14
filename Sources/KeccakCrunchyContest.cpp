/*
KeccakTools

The Keccak sponge function, designed by Guido Bertoni, Joan Daemen,
Michaël Peeters and Gilles Van Assche. For more information, feedback or
questions, please refer to our website: http://keccak.noekeon.org/

Implementation by the designers,
hereby denoted as "the implementer".

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#include <string.h>
#include "Keccak.h"
#include "KeccakCrunchyContest.h"

int verifyPreimageChallenge(int r, int c, int nr, const UINT8 image[], int start, const UINT8 preimage[], const int preimageLength)
{
    const int imageLength = 80;
    UINT8 output[10];
    ReducedRoundKeccak keccakRR(r, c, start, nr);
    keccakRR.absorb(preimage, preimageLength);
    keccakRR.squeeze(output, imageLength);
    cout << "Preimage challenge on " << keccakRR << ": ";
    if (memcmp(output, image, (imageLength+7)/8) == 0) {
        cout << "OK!" << endl;
        return 1;
    }
    else {
        cout << "failed." << endl;
        return 0;
    }
}

int verifyCollisionChallenge(int r, int c, int nr, int start, const UINT8 input1[], const int input1Length, const UINT8 input2[], const int input2Length)
{
    const int imageLength = 160;
    UINT8 output1[20], output2[20];
    {
        ReducedRoundKeccak keccakRR(r, c, start, nr);
        keccakRR.absorb(input1, input1Length);
        keccakRR.squeeze(output1, imageLength);
    }
    {
        ReducedRoundKeccak keccakRR(r, c, start, nr);
        keccakRR.absorb(input2, input2Length);
        keccakRR.squeeze(output2, imageLength);
        cout << "Collision challenge on " << keccakRR << ": ";
    }
    if (memcmp(output1, output2, (imageLength+7)/8) == 0) {
        cout << "OK!" << endl;
        return 1;
    }
    else {
        cout << "failed." << endl;
        return 0;
    }
}

void verifyPreimageChallenges()
{
    int counter = 0;
    counter += verifyPreimageChallenge(
        // Keccak[r=40, c=160, rounds=1]: preimage challenge
        40, 160, 1, (const UINT8*)"\xe9\xf5\x7f\x02\xa9\xb0\xeb\xd8\x44\x98",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=240, c=160, rounds=1]: preimage challenge
        240, 160, 1, (const UINT8*)"\xd9\xd6\xd3\xc8\x4d\x1a\xc1\xd7\x5f\x96",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=640, c=160, rounds=1]: preimage challenge
        640, 160, 1, (const UINT8*)"\x3f\x41\x9f\x88\x1c\x42\xcf\xfc\x5f\xd7",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=1440, c=160, rounds=1]: preimage challenge
        1440, 160, 1, (const UINT8*)"\x0f\x0a\xf7\x07\x4b\x6a\xbd\x48\x6f\x80",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=40, c=160, rounds=2]: preimage challenge
        40, 160, 2, (const UINT8*)"\x02\x4a\x55\x18\xe1\xe9\x5d\xb5\x32\x19",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=240, c=160, rounds=2]: preimage challenge
        240, 160, 2, (const UINT8*)"\x7a\xb8\x98\x1a\xda\x8f\xdb\x60\xae\xfd",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=640, c=160, rounds=2]: preimage challenge
        640, 160, 2, (const UINT8*)"\x82\x8d\x4d\x09\x05\x0e\x06\x35\x07\x5e",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=1440, c=160, rounds=2]: preimage challenge
        1440, 160, 2, (const UINT8*)"\x63\x90\x22\x0e\x7b\x5d\x32\x84\xd2\x3e",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=40, c=160, rounds=3]: preimage challenge
        40, 160, 3, (const UINT8*)"\xd8\xed\x85\x69\x2a\xfb\xee\x4c\x99\xce",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=240, c=160, rounds=3]: preimage challenge
        240, 160, 3, (const UINT8*)"\x5c\x9d\x5e\x4b\x38\x5e\x9c\x4f\x8e\x2e",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=640, c=160, rounds=3]: preimage challenge
        640, 160, 3, (const UINT8*)"\x00\x7b\xb5\xc5\x99\x80\x66\x0e\x02\x93",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=1440, c=160, rounds=3]: preimage challenge
        1440, 160, 3, (const UINT8*)"\x06\x25\xa3\x46\x28\xc0\xcf\xe7\x6c\x75",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=40, c=160, rounds=4]: preimage challenge
        40, 160, 4, (const UINT8*)"\x74\x2c\x7e\x3c\xd9\x46\x1d\x0d\x03\x4e",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=240, c=160, rounds=4]: preimage challenge
        240, 160, 4, (const UINT8*)"\x0d\xd2\x5e\x6d\xe2\x9a\x42\xad\xb3\x58",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=640, c=160, rounds=4]: preimage challenge
        640, 160, 4, (const UINT8*)"\x75\x1a\x16\xe5\xe4\x95\xe1\xe2\xff\x22",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=1440, c=160, rounds=4]: preimage challenge
        1440, 160, 4, (const UINT8*)"\x7d\xaa\xd8\x07\xf8\x50\x6c\x9c\x02\x76",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=40, c=160, rounds=5]: preimage challenge
        40, 160, 5, (const UINT8*)"\xe0\x53\xf9\x64\x4f\xaa\xb1\xda\x31\x1b",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=240, c=160, rounds=5]: preimage challenge
        240, 160, 5, (const UINT8*)"\x8d\xf4\x44\x09\xb4\x6f\xb8\xc6\x1b\xc4",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=640, c=160, rounds=5]: preimage challenge
        640, 160, 5, (const UINT8*)"\x6e\xf2\x61\x6f\xeb\xb9\x9b\x1f\x70\xed",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=1440, c=160, rounds=5]: preimage challenge
        1440, 160, 5, (const UINT8*)"\x65\x3b\xc0\xf8\x7d\x26\x4f\x08\x57\xd0",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=40, c=160, rounds=6]: preimage challenge
        40, 160, 6, (const UINT8*)"\xe5\x1c\x00\xc4\x8e\xd5\xdb\x07\x02\xb3",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=240, c=160, rounds=6]: preimage challenge
        240, 160, 6, (const UINT8*)"\x57\x16\xe7\x01\xef\x67\xcc\x04\x48\xb0",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=640, c=160, rounds=6]: preimage challenge
        640, 160, 6, (const UINT8*)"\x5f\x9e\x63\x88\x4f\x2e\x94\xf1\xa1\x0e",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=1440, c=160, rounds=6]: preimage challenge
        1440, 160, 6, (const UINT8*)"\xd6\x05\x33\x5e\xdc\xe7\xd2\xca\xf4\x10",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=40, c=160, rounds=7]: preimage challenge
        40, 160, 7, (const UINT8*)"\x95\x93\x25\xc5\x67\x73\xa7\x4a\x43\xc6",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=240, c=160, rounds=7]: preimage challenge
        240, 160, 7, (const UINT8*)"\x9c\xec\xce\x92\x93\x8a\xea\xba\x26\xaf",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=640, c=160, rounds=7]: preimage challenge
        640, 160, 7, (const UINT8*)"\xa4\xc1\x35\x21\x90\x12\xaa\xc8\x08\xed",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=1440, c=160, rounds=7]: preimage challenge
        1440, 160, 7, (const UINT8*)"\x5e\x0d\x17\x9c\x50\xc2\x93\x0c\x0d\x76",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=40, c=160, rounds=8]: preimage challenge
        40, 160, 8, (const UINT8*)"\x05\x4d\xda\xf1\xb9\xb5\x9b\x9a\x60\xbf",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=240, c=160, rounds=8]: preimage challenge
        240, 160, 8, (const UINT8*)"\x19\xc2\xd8\xff\x69\xe5\x66\xa5\x07\xc9",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=640, c=160, rounds=8]: preimage challenge
        640, 160, 8, (const UINT8*)"\xf4\x83\x5d\x80\x2a\xab\xc5\xbe\x75\x8e",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=1440, c=160, rounds=8]: preimage challenge
        1440, 160, 8, (const UINT8*)"\x34\xe1\x81\x23\x29\xd5\xe8\x9d\x67\x1a",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=40, c=160, rounds=9]: preimage challenge
        40, 160, 9, (const UINT8*)"\x5e\xd1\xa9\xc1\x84\xeb\x72\xb9\x45\x46",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=240, c=160, rounds=9]: preimage challenge
        240, 160, 9, (const UINT8*)"\x78\xd6\x58\xde\xc5\x01\xee\xd6\x3b\x1e",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=640, c=160, rounds=9]: preimage challenge
        640, 160, 9, (const UINT8*)"\x2e\xdd\x24\x58\x7f\x22\x5c\x69\x6e\x61",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=1440, c=160, rounds=9]: preimage challenge
        1440, 160, 9, (const UINT8*)"\xca\x18\x6a\x0f\xe1\x26\xed\xbe\x2c\xa6",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=40, c=160, rounds=10]: preimage challenge
        40, 160, 10, (const UINT8*)"\xc3\x8f\x61\x8f\x53\xa9\x6e\x4f\xfd\x53",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=240, c=160, rounds=10]: preimage challenge
        240, 160, 10, (const UINT8*)"\x46\x68\x1a\x4a\x3a\x97\x5b\x16\x2a\xc4",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=640, c=160, rounds=10]: preimage challenge
        640, 160, 10, (const UINT8*)"\xb8\x6d\xb6\x0f\xf7\x23\x18\x76\x6e\xef",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=1440, c=160, rounds=10]: preimage challenge
        1440, 160, 10, (const UINT8*)"\xdf\x7b\xf3\x01\x7c\xd3\x22\xa4\x6c\x31",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=40, c=160, rounds=11]: preimage challenge
        40, 160, 11, (const UINT8*)"\x19\xf8\xe6\xbc\x5d\x71\x41\x77\x65\x95",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=240, c=160, rounds=11]: preimage challenge
        240, 160, 11, (const UINT8*)"\x12\x9e\x94\x0f\x63\x43\x00\xf6\xb4\x14",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=640, c=160, rounds=11]: preimage challenge
        640, 160, 11, (const UINT8*)"\xa2\x49\x0a\x3e\x68\xd5\xd0\x2d\xd4\xaa",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=1440, c=160, rounds=11]: preimage challenge
        1440, 160, 11, (const UINT8*)"\x69\xc9\x4f\x0a\xe8\x30\x40\x26\xb3\xda",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=40, c=160, rounds=12]: preimage challenge
        40, 160, 12, (const UINT8*)"\x20\x68\x65\xeb\x08\xb4\x2a\x66\x63\xe1",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=240, c=160, rounds=12]: preimage challenge
        240, 160, 12, (const UINT8*)"\x85\x5a\x86\x45\x96\xc5\x1c\xaf\x7d\x3d",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=640, c=160, rounds=12]: preimage challenge
        640, 160, 12, (const UINT8*)"\x68\xed\xde\x13\xa4\x79\xe1\x47\x71\xbd",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );
    counter += verifyPreimageChallenge(
        // Keccak[r=1440, c=160, rounds=12]: preimage challenge
        1440, 160, 12, (const UINT8*)"\xbf\x8c\x82\x63\xa9\x87\x59\x5b\x21\xc0",
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24 // fill in this line
     );

    cout << dec << counter << " correct preimage challenge(s)." << endl;
}

void verifyCollisionChallenges()
{
    int counter = 0;
    counter += verifyCollisionChallenge(
        // Keccak[r=40, c=160, rounds=1]: collision challenge
        40, 160, 1,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=240, c=160, rounds=1]: collision challenge
        240, 160, 1,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=640, c=160, rounds=1]: collision challenge
        640, 160, 1,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=1440, c=160, rounds=1]: collision challenge
        1440, 160, 1,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=40, c=160, rounds=2]: collision challenge
        40, 160, 2,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=240, c=160, rounds=2]: collision challenge
        240, 160, 2,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=640, c=160, rounds=2]: collision challenge
        640, 160, 2,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=1440, c=160, rounds=2]: collision challenge
        1440, 160, 2,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=40, c=160, rounds=3]: collision challenge
        40, 160, 3,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=240, c=160, rounds=3]: collision challenge
        240, 160, 3,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=640, c=160, rounds=3]: collision challenge
        640, 160, 3,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=1440, c=160, rounds=3]: collision challenge
        1440, 160, 3,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=40, c=160, rounds=4]: collision challenge
        40, 160, 4,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=240, c=160, rounds=4]: collision challenge
        240, 160, 4,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
    );
    counter += verifyCollisionChallenge(
        // Keccak[r=640, c=160, rounds=4]: collision challenge
        640, 160, 4,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
    );
    counter += verifyCollisionChallenge(
        // Keccak[r=1440, c=160, rounds=4]: collision challenge
        1440, 160, 4,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=40, c=160, rounds=5]: collision challenge
        40, 160, 5,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=240, c=160, rounds=5]: collision challenge
        240, 160, 5,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=640, c=160, rounds=5]: collision challenge
        640, 160, 5,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=1440, c=160, rounds=5]: collision challenge
        1440, 160, 5,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=40, c=160, rounds=6]: collision challenge
        40, 160, 6,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=240, c=160, rounds=6]: collision challenge
        240, 160, 6,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=640, c=160, rounds=6]: collision challenge
        640, 160, 6,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=1440, c=160, rounds=6]: collision challenge
        1440, 160, 6,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=40, c=160, rounds=7]: collision challenge
        40, 160, 7,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=240, c=160, rounds=7]: collision challenge
        240, 160, 7,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=640, c=160, rounds=7]: collision challenge
        640, 160, 7,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=1440, c=160, rounds=7]: collision challenge
        1440, 160, 7,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=40, c=160, rounds=8]: collision challenge
        40, 160, 8,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=240, c=160, rounds=8]: collision challenge
        240, 160, 8,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=640, c=160, rounds=8]: collision challenge
        640, 160, 8,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=1440, c=160, rounds=8]: collision challenge
        1440, 160, 8,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=40, c=160, rounds=9]: collision challenge
        40, 160, 9,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=240, c=160, rounds=9]: collision challenge
        240, 160, 9,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=640, c=160, rounds=9]: collision challenge
        640, 160, 9,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=1440, c=160, rounds=9]: collision challenge
        1440, 160, 9,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=40, c=160, rounds=10]: collision challenge
        40, 160, 10,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=240, c=160, rounds=10]: collision challenge
        240, 160, 10,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=640, c=160, rounds=10]: collision challenge
        640, 160, 10,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=1440, c=160, rounds=10]: collision challenge
        1440, 160, 10,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=40, c=160, rounds=11]: collision challenge
        40, 160, 11,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=240, c=160, rounds=11]: collision challenge
        240, 160, 11,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=640, c=160, rounds=11]: collision challenge
        640, 160, 11,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=1440, c=160, rounds=11]: collision challenge
        1440, 160, 11,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=40, c=160, rounds=12]: collision challenge
        40, 160, 12,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=240, c=160, rounds=12]: collision challenge
        240, 160, 12,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=640, c=160, rounds=12]: collision challenge
        640, 160, 12,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );
    counter += verifyCollisionChallenge(
        // Keccak[r=1440, c=160, rounds=12]: collision challenge
        1440, 160, 12,
        0, // fill in this line with the start round index (0=first)
        (const UINT8*)"???", 24, // fill in this line
        (const UINT8*)"!!!!", 32 // and this line with a different input
     );

    cout << dec << counter << " correct collision challenge(s)." << endl;
}

void verifyChallenges()
{
    verifyPreimageChallenges();
    verifyCollisionChallenges();
}
