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
 * genTables.c
 *
 *  Last modified: 2021-08-01
 *  Author: Yufeng Tang, Zheng Gong
 *
 */
#include "genTables.h"

u32 TypeII[10][16][256];//Type II
u32 TypeII_R2[16][256][256];//
u32 Type_mask[16][256];
u32 TypeIII[9][16][256];//Type III
u8 TypeIV[9][96][16][16];

// static u32 nibble[16] = {0x01, 0x02, 0x0C, 0x05, 0x07, 0x08, 0x0A, 0x0F, 0x04, 0x0D, 0x0B, 0x0E, 0x09, 0x06, 0x00, 0x03};
// static u32 nibble_inv[16] = {0x0e, 0x00, 0x01, 0x0f, 0x08, 0x03, 0x0d, 0x04, 0x05, 0x0c, 0x06, 0x0a, 0x02, 0x09, 0x0b, 0x07}; 
static u8 nibble[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f};
static u8 nibble_inv[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f};


void printstate(unsigned char * in){
        for(int i = 0; i < 16; i++) {
                printf("%.2X", in[i]);

        }
        printf("\n");

        return;
}

void computeTables (u8 expandedKey[176])
{
    M8 L[9][16];
    M8 L_inv[9][16];
    M8 MaskMat[16];
    M8 MaskMat_inv[16];
    M32 MB[9][4];
    M32 MB_inv[9][4];
    for(int i = 0; i < 9; i++)
    {
        for(int j = 0; j < 16; j++)
        {
            genMatpairM8(&L[i][j], &L_inv[i][j]);
        }
    }
    for(int j = 0; j < 16; j++)
    {
        genMatpairM8(&MaskMat[j], &MaskMat_inv[j]);
    }
    for(int i = 0; i < 9; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            genMatpairM32(&MB[i][j], &MB_inv[i][j]);
        }
    }

    u32 Tyi[4][256];
    for (int x = 0; x < 256; x++)
    {
      Tyi[0][x] = (gMul(2, x) << 24) | (x << 16) | (x << 8) | gMul(3, x);
      Tyi[1][x] = (gMul(3, x) << 24) | (gMul(2, x) << 16) | (x << 8) | x;
      Tyi[2][x] = (x << 24) | (gMul(3, x) << 16) | (gMul(2, x) << 8) | x;
      Tyi[3][x] = (x << 24) | (x << 16) | (gMul(3, x) << 8) | gMul(2, x);
    }

    M32 Out_L[9][4];
    for(int i = 0; i < 9; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            MatrixcomM8to32(L[i][4 * j], L[i][4 * j + 1], L[i][4 * j + 2], L[i][4 * j + 3], &Out_L[i][j]);
        }
    }
    M32 Out_Mask[4];
    for(int j = 0; j < 4; j++)
    {
        MatrixcomM8to32(MaskMat[4 * j], MaskMat[4 * j + 1], MaskMat[4 * j + 2], MaskMat[4 * j + 3], &Out_Mask[j]);
    }
    
    int columnindex[]={0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
    //Round 1
    shiftRows (expandedKey + 16 * 0);
    printf("key: %.2x ", expandedKey[16 * 0 + 0]);
    printf("%.2x\n", expandedKey[16 * 0 + 1]);
    InitRandom(((unsigned int)time(NULL)));
    for(int j = 0; j < 16; j++)//type_II
    {
        u8 temp_u8;
        u32 temp_u32;
        u32 temp_mask;
        for(int x = 0; x < 256; x++)
        {
            temp_u8 = SBox[x ^ expandedKey[16 * 0 + j]];
            temp_u32 = Tyi[j % 4][temp_u8];
            temp_mask = random();
            temp_u32 ^= temp_mask;
            temp_u32 = MatMulNumM32(MB[0][columnindex[j]], temp_u32);
            TypeII[0][j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
        
            //table mask
            temp_u32 = MatMulNumM32(Out_Mask[columnindex[j]], temp_mask);
            Type_mask[j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
        }
    }
    for(int j = 0; j < 16; j++)//type_III
    {
        u8 temp_u8;
        u32 temp_u32;
        int shiftbit[]={24, 16, 8, 0};
        for(int x = 0; x < 256; x++)
        {
            temp_u8 = x;
            temp_u8 = (nibble_inv[(temp_u8 & 0xf0) >> 4] << 4) | (nibble_inv[(temp_u8 & 0x0f)]); 
            temp_u32 = temp_u8;
            temp_u32 = temp_u32 << shiftbit[j % 4];
            temp_u32 = MatMulNumM32(MB_inv[0][columnindex[j]], temp_u32);
            temp_u32 = MatMulNumM32(Out_L[0][columnindex[j]], temp_u32);
            TypeIII[0][j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
        }
    }

    int shiftindex[]={0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11};
    //Round 2
    shiftRows (expandedKey + 16 * 1);
    for(int j = 0; j < 16; j++)
    {
        u8 temp_u8_left;
        u8 temp_u8_right;
        u8 temp_u8;
        u32 temp_u32;
        for(int x = 0; x < 256; x++)
        {
            for(int y = 0; y < 256; y++)
            {
                temp_u8_left = x;
                temp_u8_right = y;
                temp_u8_left = (nibble_inv[(temp_u8_left & 0xf0) >> 4] << 4) | (nibble_inv[(temp_u8_left & 0x0f)]);
                temp_u8_right = (nibble_inv[(temp_u8_right & 0xf0) >> 4] << 4) | (nibble_inv[(temp_u8_right & 0x0f)]);
                
                temp_u8_left = MatMulNumM8(L_inv[0][shiftindex[j]], temp_u8_left);
                temp_u8_right = MatMulNumM8(MaskMat_inv[shiftindex[j]], temp_u8_right);
                temp_u8 = temp_u8_left ^ temp_u8_right;
                temp_u8 = SBox[temp_u8 ^ expandedKey[16 * 1 + j]];
                temp_u32 = Tyi[j % 4][temp_u8];
                temp_u32 = MatMulNumM32(MB[1][columnindex[j]], temp_u32);
                TypeII_R2[j][x][y] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
            }
        }
    }
    for(int j = 0; j < 16; j++)//type_III
    {
        u8 temp_u8;
        u32 temp_u32;
        int shiftbit[]={24, 16, 8, 0};
        for(int x = 0; x < 256; x++)
        {
            temp_u8 = x;
            temp_u8 = (nibble_inv[(temp_u8 & 0xf0) >> 4] << 4) | (nibble_inv[(temp_u8 & 0x0f)]);
            temp_u32 = temp_u8;
            temp_u32 = temp_u32 << shiftbit[j % 4];
            temp_u32 = MatMulNumM32(MB_inv[1][columnindex[j]], temp_u32);
            temp_u32 = MatMulNumM32(Out_L[1][columnindex[j]], temp_u32);
            TypeIII[1][j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
        }
    }

    //Round 3-9
    for (int i = 2; i < 9; i++)//Type_II
    {
        shiftRows (expandedKey + 16 * i);
        for(int j = 0; j < 16; j++)
        {
            u8 temp_u8;
            u32 temp_u32;
            for(int x = 0; x < 256; x++)
            {
                temp_u8 = x;
                temp_u8 = (nibble_inv[(temp_u8 & 0xf0) >> 4] << 4) | (nibble_inv[(temp_u8 & 0x0f)]);
                temp_u8 = MatMulNumM8(L_inv[i - 1][shiftindex[j]], temp_u8);
                temp_u8 = SBox[temp_u8 ^ expandedKey[16 * i + j]];
                temp_u32 = Tyi[j % 4][temp_u8];
                temp_u32 = MatMulNumM32(MB[i][columnindex[j]], temp_u32);
                TypeII[i][j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
            }
        }
    
        for(int j = 0; j < 16; j++)//type_III
        {
            u8 temp_u8;
            u32 temp_u32;
            int shiftbit[]={24, 16, 8, 0};
            for(int x = 0; x < 256; x++)
            {
                temp_u8 = x;
                temp_u8 = (nibble_inv[(temp_u8 & 0xf0) >> 4] << 4) | (nibble_inv[(temp_u8 & 0x0f)]);
                temp_u32 = temp_u8;
                temp_u32 = temp_u32 << shiftbit[j % 4];
                temp_u32 = MatMulNumM32(MB_inv[i][columnindex[j]], temp_u32);
                temp_u32 = MatMulNumM32(Out_L[i][columnindex[j]], temp_u32);
                TypeIII[i][j][x] = (nibble[(temp_u32 & 0xf0000000) >> 28] << 28) | (nibble[(temp_u32 & 0x0f000000) >> 24] << 24) | (nibble[(temp_u32 & 0x00f00000) >> 20] << 20) | (nibble[(temp_u32 & 0x000f0000) >> 16] << 16) | (nibble[(temp_u32 & 0x0000f000) >> 12] << 12) | (nibble[(temp_u32 & 0x00000f00) >> 8] << 8) | (nibble[(temp_u32 & 0x000000f0) >> 4] << 4) | (nibble[(temp_u32 & 0x0000000f)]);
            }
        }
    }

    //Round 10
    shiftRows (expandedKey + 16 * 9);
    for(int j = 0; j < 16; j++)//type_II
    {
        u8 temp_u8;
        for(int x = 0; x < 256; x++)
        {
            temp_u8 = x;
            temp_u8 = (nibble_inv[(temp_u8 & 0xf0) >> 4] << 4) | (nibble_inv[(temp_u8 & 0x0f)]);
            temp_u8 = MatMulNumM8(L_inv[8][shiftindex[j]], temp_u8);
            temp_u8 = SBox[temp_u8 ^ expandedKey[16 * 9 + j]];
            TypeII[9][j][x] = temp_u8 ^ expandedKey[16 * 10 + j];
        }
    }

    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 96; j++)
        {
            for (int x = 0; x < 16; x++)
            {
                for (int y = 0; y < 16; y++)
                {
                    TypeIV[i][j][x][y] = nibble[nibble_inv[x] ^ nibble_inv[y]];
                }
            }
        }
    }
}
int aes_128_table_encrypt_collect (int collect[256][Trace]) 
{
  u32 a, b, c, d, aa, bb, cc, dd;
  u8 mask[16];
  u8 mask_input[16];
    //Round 1
    u8 input[16];
    int count = 0;
    for(int left = 0; left < 256; left++)
    {
        for(int right = 0; right < 256; right++)
        {
            mask_input[0] = left;
            mask_input[1] = right;
            for(int i = 2; i < 16; i++)
            {
                mask_input[i] = 0;
            }
            for (int j = 0; j < 4; j++)
            {
                a = TypeII[0][4*j + 0][input[4*j + 0]];
                b = TypeII[0][4*j + 1][input[4*j + 1]];
                c = TypeII[0][4*j + 2][input[4*j + 2]];
                d = TypeII[0][4*j + 3][input[4*j + 3]];

                aa = TypeIV[0][24*j + 0][(a >> 28) & 0xf][(b >> 28) & 0xf];
                bb = TypeIV[0][24*j + 1][(c >> 28) & 0xf][(d >> 28) & 0xf];
                cc = TypeIV[0][24*j + 2][(a >> 24) & 0xf][(b >> 24) & 0xf];
                dd = TypeIV[0][24*j + 3][(c >> 24) & 0xf][(d >> 24) & 0xf];
                input[4*j + 0] = (TypeIV[0][24*j + 4][aa][bb] << 4) | TypeIV[0][24*j + 5][cc][dd];

                aa = TypeIV[0][24*j + 6][(a >> 20) & 0xf][(b >> 20) & 0xf];
                bb = TypeIV[0][24*j + 7][(c >> 20) & 0xf][(d >> 20) & 0xf];
                cc = TypeIV[0][24*j + 8][(a >> 16) & 0xf][(b >> 16) & 0xf];
                dd = TypeIV[0][24*j + 9][(c >> 16) & 0xf][(d >> 16) & 0xf];
                input[4*j + 1] = (TypeIV[0][24*j + 10][aa][bb] << 4) | TypeIV[0][24*j + 11][cc][dd];

                aa = TypeIV[0][24*j + 12][(a >> 12) & 0xf][(b >> 12) & 0xf];
                bb = TypeIV[0][24*j + 13][(c >> 12) & 0xf][(d >> 12) & 0xf];
                cc = TypeIV[0][24*j + 14][(a >>  8) & 0xf][(b >>  8) & 0xf];
                dd = TypeIV[0][24*j + 15][(c >>  8) & 0xf][(d >>  8) & 0xf];
                input[4*j + 2] = (TypeIV[0][24*j + 16][aa][bb] << 4) | TypeIV[0][24*j + 17][cc][dd];

                aa = TypeIV[0][24*j + 18][(a >>  4) & 0xf][(b >>  4) & 0xf];
                bb = TypeIV[0][24*j + 19][(c >>  4) & 0xf][(d >>  4) & 0xf];
                cc = TypeIV[0][24*j + 20][(a >>  0) & 0xf][(b >>  0) & 0xf];
                dd = TypeIV[0][24*j + 21][(c >>  0) & 0xf][(d >>  0) & 0xf];
                input[4*j + 3] = (TypeIV[0][24*j + 22][aa][bb] << 4) | TypeIV[0][24*j + 23][cc][dd];


                a = TypeIII[0][4*j + 0][input[4*j + 0]];
                b = TypeIII[0][4*j + 1][input[4*j + 1]];
                c = TypeIII[0][4*j + 2][input[4*j + 2]];
                d = TypeIII[0][4*j + 3][input[4*j + 3]];

                aa = TypeIV[0][24*j + 0][(a >> 28) & 0xf][(b >> 28) & 0xf];
                bb = TypeIV[0][24*j + 1][(c >> 28) & 0xf][(d >> 28) & 0xf];
                cc = TypeIV[0][24*j + 2][(a >> 24) & 0xf][(b >> 24) & 0xf];
                dd = TypeIV[0][24*j + 3][(c >> 24) & 0xf][(d >> 24) & 0xf];
                input[4*j + 0] = (TypeIV[0][24*j + 4][aa][bb] << 4) | TypeIV[0][24*j + 5][cc][dd];

                aa = TypeIV[0][24*j + 6][(a >> 20) & 0xf][(b >> 20) & 0xf];
                bb = TypeIV[0][24*j + 7][(c >> 20) & 0xf][(d >> 20) & 0xf];
                cc = TypeIV[0][24*j + 8][(a >> 16) & 0xf][(b >> 16) & 0xf];
                dd = TypeIV[0][24*j + 9][(c >> 16) & 0xf][(d >> 16) & 0xf];
                input[4*j + 1] = (TypeIV[0][24*j + 10][aa][bb] << 4) | TypeIV[0][24*j + 11][cc][dd];

                aa = TypeIV[0][24*j + 12][(a >> 12) & 0xf][(b >> 12) & 0xf];
                bb = TypeIV[0][24*j + 13][(c >> 12) & 0xf][(d >> 12) & 0xf];
                cc = TypeIV[0][24*j + 14][(a >>  8) & 0xf][(b >>  8) & 0xf];
                dd = TypeIV[0][24*j + 15][(c >>  8) & 0xf][(d >>  8) & 0xf];
                input[4*j + 2] = (TypeIV[0][24*j + 16][aa][bb] << 4) | TypeIV[0][24*j + 17][cc][dd];

                aa = TypeIV[0][24*j + 18][(a >>  4) & 0xf][(b >>  4) & 0xf];
                bb = TypeIV[0][24*j + 19][(c >>  4) & 0xf][(d >>  4) & 0xf];
                cc = TypeIV[0][24*j + 20][(a >>  0) & 0xf][(b >>  0) & 0xf];
                dd = TypeIV[0][24*j + 21][(c >>  0) & 0xf][(d >>  0) & 0xf];
                input[4*j + 3] = (TypeIV[0][24*j + 22][aa][bb] << 4) | TypeIV[0][24*j + 23][cc][dd];     

                //mask
                a = Type_mask[4*j + 0][mask_input[4*j + 0]];
                b = Type_mask[4*j + 1][mask_input[4*j + 1]];
                c = Type_mask[4*j + 2][mask_input[4*j + 2]];
                d = Type_mask[4*j + 3][mask_input[4*j + 3]];

                aa = TypeIV[0][24*j + 0][(a >> 28) & 0xf][(b >> 28) & 0xf];
                bb = TypeIV[0][24*j + 1][(c >> 28) & 0xf][(d >> 28) & 0xf];
                cc = TypeIV[0][24*j + 2][(a >> 24) & 0xf][(b >> 24) & 0xf];
                dd = TypeIV[0][24*j + 3][(c >> 24) & 0xf][(d >> 24) & 0xf];
                mask[4*j + 0] = (TypeIV[0][24*j + 4][aa][bb] << 4) | TypeIV[0][24*j + 5][cc][dd];

                aa = TypeIV[0][24*j + 6][(a >> 20) & 0xf][(b >> 20) & 0xf];
                bb = TypeIV[0][24*j + 7][(c >> 20) & 0xf][(d >> 20) & 0xf];
                cc = TypeIV[0][24*j + 8][(a >> 16) & 0xf][(b >> 16) & 0xf];
                dd = TypeIV[0][24*j + 9][(c >> 16) & 0xf][(d >> 16) & 0xf];
                mask[4*j + 1] = (TypeIV[0][24*j + 10][aa][bb] << 4) | TypeIV[0][24*j + 11][cc][dd];

                aa = TypeIV[0][24*j + 12][(a >> 12) & 0xf][(b >> 12) & 0xf];
                bb = TypeIV[0][24*j + 13][(c >> 12) & 0xf][(d >> 12) & 0xf];
                cc = TypeIV[0][24*j + 14][(a >>  8) & 0xf][(b >>  8) & 0xf];
                dd = TypeIV[0][24*j + 15][(c >>  8) & 0xf][(d >>  8) & 0xf];
                mask[4*j + 2] = (TypeIV[0][24*j + 16][aa][bb] << 4) | TypeIV[0][24*j + 17][cc][dd];

                aa = TypeIV[0][24*j + 18][(a >>  4) & 0xf][(b >>  4) & 0xf];
                bb = TypeIV[0][24*j + 19][(c >>  4) & 0xf][(d >>  4) & 0xf];
                cc = TypeIV[0][24*j + 20][(a >>  0) & 0xf][(b >>  0) & 0xf];
                dd = TypeIV[0][24*j + 21][(c >>  0) & 0xf][(d >>  0) & 0xf];
                mask[4*j + 3] = (TypeIV[0][24*j + 22][aa][bb] << 4) | TypeIV[0][24*j + 23][cc][dd];
            }
            //count++;
            for(int k = 0; k < Trace; k++)
            {
                if(collect[mask[0]][k] == -1) 
                {
                    uint16_t ttval = mask_input[0];
                    ttval = ttval << 8;
                    ttval ^= mask_input[1];
                    collect[mask[0]][k] = ttval;
                    if(k == Trace - 1) return 0;
                    break;
                }
            }
        }
    }
}
void aes_128_table_encrypt (uint16_t collect_input[Trace], int middlebitstate[Trace][8]) {
  u32 a, b, c, d, aa, bb, cc, dd;
  u8 mask[16];
  u8 mask_input[16];
  u8 input[16];
    //Round 1
    for(int time = 0; time < Trace; time++)
    {  
        uint16_t ttval = collect_input[time];
        input[1] = ttval;
        ttval = ttval >> 8;
        input[0] = ttval;
        for(int i = 2; i < 16; i++)
        {
            input[i] = 0;
        }
        for (int j = 0; j < 4; j++)
        {
            a = TypeII[0][4*j + 0][input[4*j + 0]];
            b = TypeII[0][4*j + 1][input[4*j + 1]];
            c = TypeII[0][4*j + 2][input[4*j + 2]];
            d = TypeII[0][4*j + 3][input[4*j + 3]];

            aa = TypeIV[0][24*j + 0][(a >> 28) & 0xf][(b >> 28) & 0xf];
            bb = TypeIV[0][24*j + 1][(c >> 28) & 0xf][(d >> 28) & 0xf];
            cc = TypeIV[0][24*j + 2][(a >> 24) & 0xf][(b >> 24) & 0xf];
            dd = TypeIV[0][24*j + 3][(c >> 24) & 0xf][(d >> 24) & 0xf];
            input[4*j + 0] = (TypeIV[0][24*j + 4][aa][bb] << 4) | TypeIV[0][24*j + 5][cc][dd];

            aa = TypeIV[0][24*j + 6][(a >> 20) & 0xf][(b >> 20) & 0xf];
            bb = TypeIV[0][24*j + 7][(c >> 20) & 0xf][(d >> 20) & 0xf];
            cc = TypeIV[0][24*j + 8][(a >> 16) & 0xf][(b >> 16) & 0xf];
            dd = TypeIV[0][24*j + 9][(c >> 16) & 0xf][(d >> 16) & 0xf];
            input[4*j + 1] = (TypeIV[0][24*j + 10][aa][bb] << 4) | TypeIV[0][24*j + 11][cc][dd];

            aa = TypeIV[0][24*j + 12][(a >> 12) & 0xf][(b >> 12) & 0xf];
            bb = TypeIV[0][24*j + 13][(c >> 12) & 0xf][(d >> 12) & 0xf];
            cc = TypeIV[0][24*j + 14][(a >>  8) & 0xf][(b >>  8) & 0xf];
            dd = TypeIV[0][24*j + 15][(c >>  8) & 0xf][(d >>  8) & 0xf];
            input[4*j + 2] = (TypeIV[0][24*j + 16][aa][bb] << 4) | TypeIV[0][24*j + 17][cc][dd];

            aa = TypeIV[0][24*j + 18][(a >>  4) & 0xf][(b >>  4) & 0xf];
            bb = TypeIV[0][24*j + 19][(c >>  4) & 0xf][(d >>  4) & 0xf];
            cc = TypeIV[0][24*j + 20][(a >>  0) & 0xf][(b >>  0) & 0xf];
            dd = TypeIV[0][24*j + 21][(c >>  0) & 0xf][(d >>  0) & 0xf];
            input[4*j + 3] = (TypeIV[0][24*j + 22][aa][bb] << 4) | TypeIV[0][24*j + 23][cc][dd];


            a = TypeIII[0][4*j + 0][input[4*j + 0]];
            b = TypeIII[0][4*j + 1][input[4*j + 1]];
            c = TypeIII[0][4*j + 2][input[4*j + 2]];
            d = TypeIII[0][4*j + 3][input[4*j + 3]];

            aa = TypeIV[0][24*j + 0][(a >> 28) & 0xf][(b >> 28) & 0xf];
            bb = TypeIV[0][24*j + 1][(c >> 28) & 0xf][(d >> 28) & 0xf];
            cc = TypeIV[0][24*j + 2][(a >> 24) & 0xf][(b >> 24) & 0xf];
            dd = TypeIV[0][24*j + 3][(c >> 24) & 0xf][(d >> 24) & 0xf];
            input[4*j + 0] = (TypeIV[0][24*j + 4][aa][bb] << 4) | TypeIV[0][24*j + 5][cc][dd];

            aa = TypeIV[0][24*j + 6][(a >> 20) & 0xf][(b >> 20) & 0xf];
            bb = TypeIV[0][24*j + 7][(c >> 20) & 0xf][(d >> 20) & 0xf];
            cc = TypeIV[0][24*j + 8][(a >> 16) & 0xf][(b >> 16) & 0xf];
            dd = TypeIV[0][24*j + 9][(c >> 16) & 0xf][(d >> 16) & 0xf];
            input[4*j + 1] = (TypeIV[0][24*j + 10][aa][bb] << 4) | TypeIV[0][24*j + 11][cc][dd];

            aa = TypeIV[0][24*j + 12][(a >> 12) & 0xf][(b >> 12) & 0xf];
            bb = TypeIV[0][24*j + 13][(c >> 12) & 0xf][(d >> 12) & 0xf];
            cc = TypeIV[0][24*j + 14][(a >>  8) & 0xf][(b >>  8) & 0xf];
            dd = TypeIV[0][24*j + 15][(c >>  8) & 0xf][(d >>  8) & 0xf];
            input[4*j + 2] = (TypeIV[0][24*j + 16][aa][bb] << 4) | TypeIV[0][24*j + 17][cc][dd];

            aa = TypeIV[0][24*j + 18][(a >>  4) & 0xf][(b >>  4) & 0xf];
            bb = TypeIV[0][24*j + 19][(c >>  4) & 0xf][(d >>  4) & 0xf];
            cc = TypeIV[0][24*j + 20][(a >>  0) & 0xf][(b >>  0) & 0xf];
            dd = TypeIV[0][24*j + 21][(c >>  0) & 0xf][(d >>  0) & 0xf];
            input[4*j + 3] = (TypeIV[0][24*j + 22][aa][bb] << 4) | TypeIV[0][24*j + 23][cc][dd];     

            // mask
            a = Type_mask[4*j + 0][mask_input[4*j + 0]];
            b = Type_mask[4*j + 1][mask_input[4*j + 1]];
            c = Type_mask[4*j + 2][mask_input[4*j + 2]];
            d = Type_mask[4*j + 3][mask_input[4*j + 3]];

            aa = TypeIV[0][24*j + 0][(a >> 28) & 0xf][(b >> 28) & 0xf];
            bb = TypeIV[0][24*j + 1][(c >> 28) & 0xf][(d >> 28) & 0xf];
            cc = TypeIV[0][24*j + 2][(a >> 24) & 0xf][(b >> 24) & 0xf];
            dd = TypeIV[0][24*j + 3][(c >> 24) & 0xf][(d >> 24) & 0xf];
            mask[4*j + 0] = (TypeIV[0][24*j + 4][aa][bb] << 4) | TypeIV[0][24*j + 5][cc][dd];

            aa = TypeIV[0][24*j + 6][(a >> 20) & 0xf][(b >> 20) & 0xf];
            bb = TypeIV[0][24*j + 7][(c >> 20) & 0xf][(d >> 20) & 0xf];
            cc = TypeIV[0][24*j + 8][(a >> 16) & 0xf][(b >> 16) & 0xf];
            dd = TypeIV[0][24*j + 9][(c >> 16) & 0xf][(d >> 16) & 0xf];
            mask[4*j + 1] = (TypeIV[0][24*j + 10][aa][bb] << 4) | TypeIV[0][24*j + 11][cc][dd];

            aa = TypeIV[0][24*j + 12][(a >> 12) & 0xf][(b >> 12) & 0xf];
            bb = TypeIV[0][24*j + 13][(c >> 12) & 0xf][(d >> 12) & 0xf];
            cc = TypeIV[0][24*j + 14][(a >>  8) & 0xf][(b >>  8) & 0xf];
            dd = TypeIV[0][24*j + 15][(c >>  8) & 0xf][(d >>  8) & 0xf];
            mask[4*j + 2] = (TypeIV[0][24*j + 16][aa][bb] << 4) | TypeIV[0][24*j + 17][cc][dd];

            aa = TypeIV[0][24*j + 18][(a >>  4) & 0xf][(b >>  4) & 0xf];
            bb = TypeIV[0][24*j + 19][(c >>  4) & 0xf][(d >>  4) & 0xf];
            cc = TypeIV[0][24*j + 20][(a >>  0) & 0xf][(b >>  0) & 0xf];
            dd = TypeIV[0][24*j + 21][(c >>  0) & 0xf][(d >>  0) & 0xf];
            mask[4*j + 3] = (TypeIV[0][24*j + 22][aa][bb] << 4) | TypeIV[0][24*j + 23][cc][dd];
        }
        u8 result = input[0];
        for(int k = 0; k < 8; k++)
        {
            if(result & 0x80) middlebitstate[time][k] = 1;
            else middlebitstate[time][k] = 0;
            result = result << 1;
        }
    }

/*
    //Round 2
    shiftRows (input);
    shiftRows (mask);
    for (int j = 0; j < 4; j++)
    {
      int i = 1;
      a = TypeII_R2[4*j + 0][input[4*j + 0]][mask[4*j + 0]];
      b = TypeII_R2[4*j + 1][input[4*j + 1]][mask[4*j + 1]];
      c = TypeII_R2[4*j + 2][input[4*j + 2]][mask[4*j + 2]];
      d = TypeII_R2[4*j + 3][input[4*j + 3]][mask[4*j + 3]];

      aa = TypeIV[i][24*j + 0][(a >> 28) & 0xf][(b >> 28) & 0xf];
      bb = TypeIV[i][24*j + 1][(c >> 28) & 0xf][(d >> 28) & 0xf];
      cc = TypeIV[i][24*j + 2][(a >> 24) & 0xf][(b >> 24) & 0xf];
      dd = TypeIV[i][24*j + 3][(c >> 24) & 0xf][(d >> 24) & 0xf];
      input[4*j + 0] = (TypeIV[i][24*j + 4][aa][bb] << 4) | TypeIV[i][24*j + 5][cc][dd];

      aa = TypeIV[i][24*j + 6][(a >> 20) & 0xf][(b >> 20) & 0xf];
      bb = TypeIV[i][24*j + 7][(c >> 20) & 0xf][(d >> 20) & 0xf];
      cc = TypeIV[i][24*j + 8][(a >> 16) & 0xf][(b >> 16) & 0xf];
      dd = TypeIV[i][24*j + 9][(c >> 16) & 0xf][(d >> 16) & 0xf];
      input[4*j + 1] = (TypeIV[i][24*j + 10][aa][bb] << 4) | TypeIV[i][24*j + 11][cc][dd];

      aa = TypeIV[i][24*j + 12][(a >> 12) & 0xf][(b >> 12) & 0xf];
      bb = TypeIV[i][24*j + 13][(c >> 12) & 0xf][(d >> 12) & 0xf];
      cc = TypeIV[i][24*j + 14][(a >>  8) & 0xf][(b >>  8) & 0xf];
      dd = TypeIV[i][24*j + 15][(c >>  8) & 0xf][(d >>  8) & 0xf];
      input[4*j + 2] = (TypeIV[i][24*j + 16][aa][bb] << 4) | TypeIV[i][24*j + 17][cc][dd];

      aa = TypeIV[i][24*j + 18][(a >>  4) & 0xf][(b >>  4) & 0xf];
      bb = TypeIV[i][24*j + 19][(c >>  4) & 0xf][(d >>  4) & 0xf];
      cc = TypeIV[i][24*j + 20][(a >>  0) & 0xf][(b >>  0) & 0xf];
      dd = TypeIV[i][24*j + 21][(c >>  0) & 0xf][(d >>  0) & 0xf];
      input[4*j + 3] = (TypeIV[i][24*j + 22][aa][bb] << 4) | TypeIV[i][24*j + 23][cc][dd];


      a = TypeIII[i][4*j + 0][input[4*j + 0]];
      b = TypeIII[i][4*j + 1][input[4*j + 1]];
      c = TypeIII[i][4*j + 2][input[4*j + 2]];
      d = TypeIII[i][4*j + 3][input[4*j + 3]];

      aa = TypeIV[i][24*j + 0][(a >> 28) & 0xf][(b >> 28) & 0xf];
      bb = TypeIV[i][24*j + 1][(c >> 28) & 0xf][(d >> 28) & 0xf];
      cc = TypeIV[i][24*j + 2][(a >> 24) & 0xf][(b >> 24) & 0xf];
      dd = TypeIV[i][24*j + 3][(c >> 24) & 0xf][(d >> 24) & 0xf];
      input[4*j + 0] = (TypeIV[i][24*j + 4][aa][bb] << 4) | TypeIV[i][24*j + 5][cc][dd];

      aa = TypeIV[i][24*j + 6][(a >> 20) & 0xf][(b >> 20) & 0xf];
      bb = TypeIV[i][24*j + 7][(c >> 20) & 0xf][(d >> 20) & 0xf];
      cc = TypeIV[i][24*j + 8][(a >> 16) & 0xf][(b >> 16) & 0xf];
      dd = TypeIV[i][24*j + 9][(c >> 16) & 0xf][(d >> 16) & 0xf];
      input[4*j + 1] = (TypeIV[i][24*j + 10][aa][bb] << 4) | TypeIV[i][24*j + 11][cc][dd];

      aa = TypeIV[i][24*j + 12][(a >> 12) & 0xf][(b >> 12) & 0xf];
      bb = TypeIV[i][24*j + 13][(c >> 12) & 0xf][(d >> 12) & 0xf];
      cc = TypeIV[i][24*j + 14][(a >>  8) & 0xf][(b >>  8) & 0xf];
      dd = TypeIV[i][24*j + 15][(c >>  8) & 0xf][(d >>  8) & 0xf];
      input[4*j + 2] = (TypeIV[i][24*j + 16][aa][bb] << 4) | TypeIV[i][24*j + 17][cc][dd];

      aa = TypeIV[i][24*j + 18][(a >>  4) & 0xf][(b >>  4) & 0xf];
      bb = TypeIV[i][24*j + 19][(c >>  4) & 0xf][(d >>  4) & 0xf];
      cc = TypeIV[i][24*j + 20][(a >>  0) & 0xf][(b >>  0) & 0xf];
      dd = TypeIV[i][24*j + 21][(c >>  0) & 0xf][(d >>  0) & 0xf];
      input[4*j + 3] = (TypeIV[i][24*j + 22][aa][bb] << 4) | TypeIV[i][24*j + 23][cc][dd];
      
    }
  
  for (int i = 2; i < 9; i++) {
    shiftRows (input);
    for (int j = 0; j < 4; j++)
    {
      a = TypeII[i][4*j + 0][input[4*j + 0]];
      b = TypeII[i][4*j + 1][input[4*j + 1]];
      c = TypeII[i][4*j + 2][input[4*j + 2]];
      d = TypeII[i][4*j + 3][input[4*j + 3]];

      aa = TypeIV[i][24*j + 0][(a >> 28) & 0xf][(b >> 28) & 0xf];
      bb = TypeIV[i][24*j + 1][(c >> 28) & 0xf][(d >> 28) & 0xf];
      cc = TypeIV[i][24*j + 2][(a >> 24) & 0xf][(b >> 24) & 0xf];
      dd = TypeIV[i][24*j + 3][(c >> 24) & 0xf][(d >> 24) & 0xf];
      input[4*j + 0] = (TypeIV[i][24*j + 4][aa][bb] << 4) | TypeIV[i][24*j + 5][cc][dd];

      aa = TypeIV[i][24*j + 6][(a >> 20) & 0xf][(b >> 20) & 0xf];
      bb = TypeIV[i][24*j + 7][(c >> 20) & 0xf][(d >> 20) & 0xf];
      cc = TypeIV[i][24*j + 8][(a >> 16) & 0xf][(b >> 16) & 0xf];
      dd = TypeIV[i][24*j + 9][(c >> 16) & 0xf][(d >> 16) & 0xf];
      input[4*j + 1] = (TypeIV[i][24*j + 10][aa][bb] << 4) | TypeIV[i][24*j + 11][cc][dd];

      aa = TypeIV[i][24*j + 12][(a >> 12) & 0xf][(b >> 12) & 0xf];
      bb = TypeIV[i][24*j + 13][(c >> 12) & 0xf][(d >> 12) & 0xf];
      cc = TypeIV[i][24*j + 14][(a >>  8) & 0xf][(b >>  8) & 0xf];
      dd = TypeIV[i][24*j + 15][(c >>  8) & 0xf][(d >>  8) & 0xf];
      input[4*j + 2] = (TypeIV[i][24*j + 16][aa][bb] << 4) | TypeIV[i][24*j + 17][cc][dd];

      aa = TypeIV[i][24*j + 18][(a >>  4) & 0xf][(b >>  4) & 0xf];
      bb = TypeIV[i][24*j + 19][(c >>  4) & 0xf][(d >>  4) & 0xf];
      cc = TypeIV[i][24*j + 20][(a >>  0) & 0xf][(b >>  0) & 0xf];
      dd = TypeIV[i][24*j + 21][(c >>  0) & 0xf][(d >>  0) & 0xf];
      input[4*j + 3] = (TypeIV[i][24*j + 22][aa][bb] << 4) | TypeIV[i][24*j + 23][cc][dd];


      a = TypeIII[i][4*j + 0][input[4*j + 0]];
      b = TypeIII[i][4*j + 1][input[4*j + 1]];
      c = TypeIII[i][4*j + 2][input[4*j + 2]];
      d = TypeIII[i][4*j + 3][input[4*j + 3]];

      aa = TypeIV[i][24*j + 0][(a >> 28) & 0xf][(b >> 28) & 0xf];
      bb = TypeIV[i][24*j + 1][(c >> 28) & 0xf][(d >> 28) & 0xf];
      cc = TypeIV[i][24*j + 2][(a >> 24) & 0xf][(b >> 24) & 0xf];
      dd = TypeIV[i][24*j + 3][(c >> 24) & 0xf][(d >> 24) & 0xf];
      input[4*j + 0] = (TypeIV[i][24*j + 4][aa][bb] << 4) | TypeIV[i][24*j + 5][cc][dd];

      aa = TypeIV[i][24*j + 6][(a >> 20) & 0xf][(b >> 20) & 0xf];
      bb = TypeIV[i][24*j + 7][(c >> 20) & 0xf][(d >> 20) & 0xf];
      cc = TypeIV[i][24*j + 8][(a >> 16) & 0xf][(b >> 16) & 0xf];
      dd = TypeIV[i][24*j + 9][(c >> 16) & 0xf][(d >> 16) & 0xf];
      input[4*j + 1] = (TypeIV[i][24*j + 10][aa][bb] << 4) | TypeIV[i][24*j + 11][cc][dd];

      aa = TypeIV[i][24*j + 12][(a >> 12) & 0xf][(b >> 12) & 0xf];
      bb = TypeIV[i][24*j + 13][(c >> 12) & 0xf][(d >> 12) & 0xf];
      cc = TypeIV[i][24*j + 14][(a >>  8) & 0xf][(b >>  8) & 0xf];
      dd = TypeIV[i][24*j + 15][(c >>  8) & 0xf][(d >>  8) & 0xf];
      input[4*j + 2] = (TypeIV[i][24*j + 16][aa][bb] << 4) | TypeIV[i][24*j + 17][cc][dd];

      aa = TypeIV[i][24*j + 18][(a >>  4) & 0xf][(b >>  4) & 0xf];
      bb = TypeIV[i][24*j + 19][(c >>  4) & 0xf][(d >>  4) & 0xf];
      cc = TypeIV[i][24*j + 20][(a >>  0) & 0xf][(b >>  0) & 0xf];
      dd = TypeIV[i][24*j + 21][(c >>  0) & 0xf][(d >>  0) & 0xf];
      input[4*j + 3] = (TypeIV[i][24*j + 22][aa][bb] << 4) | TypeIV[i][24*j + 23][cc][dd];
      
    }
  }
    //Round 10
    shiftRows(input);
    for (int j = 0; j < 16; j++) {
        input[j] = TypeII[9][j][input[j]];
    }

    for (int i = 0; i < 16; i++)
        output[i] = input[i];
*/
}