/*
 * Copyright (c) 2021, South China Normal University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in
 *   the documentation and/or other materials provided with the
 *   distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
/*
 * main.c
 *
 *  Last modified: 2021-08-01
 *  Author: Yufeng Tang, Zheng Gong
 *
 */
#include "genTables.h"

int main(int argc, char * argv[])
{
    #define bit 8
    u8 expandedKey[176];
    u8 key[16] = {0x7b, 0, 0, 0, 0, 0xc3};
    uint64_t begin;
    uint64_t end;
    uint64_t ans = 0;
    expandKey (key, expandedKey);
    computeTables(expandedKey);    

    int correct_times = 0;
    int collect[256][Trace];
    int middlebitstate[Trace][8] = {0};
    
    memset(collect, -1, sizeof(collect));
    aes_128_table_encrypt_collect(collect);

    int index = 0;
    for(int i = 0; i < 256; i++)
    {
        if(collect[i][Trace - 1] != -1) 
        {
            index = i;
            break;
        }
    }
    uint16_t plain[Trace] = {0};
    for(int p = 0; p < Trace; p++)
    {
        plain[p] = collect[index][p];
    }
    uint8_t plaintext[Trace][2] = {0};
    for(int i = 0; i < Trace; i++)
    {
        uint16_t ttval = plain[i];
        plaintext[i][1] =  ttval;
        ttval = ttval >> 8;
        plaintext[i][0] =  ttval;
    }
    aes_128_table_encrypt(plain, middlebitstate);

    uint8_t identM8[8]={0x80,0x40,0x20,0x10,0x08,0x04,0x02,0x01};
    u8 k1 = 0, k2 = 0;
    double keystate[8][9]={0.0};
    double maxkeystate[9]={0.0};
    u8 keyguess[2] = {0};
    uint8_t in1 = 0, in2 = 0;
    uint8_t cipher = 0;
    for(int left = 0; left<256; left++)//key
    {
        k1 = left;
        for(int right = 0; right<256; right++)//key
        {
            k2 = right;
            for(int r=0;r<8;r++)//bit
            {
                double split0[bit]={0.0};
                double split1[bit]={0.0};
                int count0=0;
                int count1=0;
                for(int i=0; i<Trace; i++)//plaintext
                {
                    in1 = plaintext[i][0] ^ k1;
                    in2 = plaintext[i][1] ^ k2;
                    in1 = SBox[in1];
                    in2 = SBox[in2];
                    cipher = gMul(2, in1) ^ gMul(3, in2);
                    if(cipher & identM8[r])
                    {
                        for(int g=0;g<bit;g++)
                        {
                            split1[g]+=middlebitstate[i][g];
                        }
                        count1++;
                    }
                    else
                    {
                        for(int g=0;g<bit;g++)
                        {
                            split0[g]+=middlebitstate[i][g];
                        }
                        count0++;
                    }
                }
                for(int y=0;y<bit;y++)
                {
                    if(count0 != 0) split0[y]=split0[y]/count0;
                    if(count1 != 0) split1[y]=split1[y]/count1;
                }
                for(int y=0;y<bit;y++)
                {
                    if(split0[y]>=split1[y]) keystate[r][y]=split0[y]-split1[y];
                    else keystate[r][y]=split1[y]-split0[y];
                }
            }
            int maxord=0;
            double maxvalue=0.0;
            for(int a=0; a<8; a++)
            {
                double size=0.0;
                for(int b=0; b<bit; b++)
                {
                    if (keystate[a][b]>size) size=keystate[a][b];
                }
                keystate[a][bit]=size;
            }
            for(int a=0; a<8; a++)
            {
                if(keystate[a][bit] > maxvalue) 
                {
                    maxord = a;
                    maxvalue = keystate[a][bit];
                }
            }
            if(maxvalue > maxkeystate[bit])
            {    
                for(int a=0; a<9; a++)
                {
                    maxkeystate[a] = keystate[maxord][a];
                    keyguess[0] = k1;
                    keyguess[1] = k2;
                }
            }
        }
    }
    printf("key guess: %.2x %.2x \n", keyguess[0], keyguess[1]);
    return 0;
}
