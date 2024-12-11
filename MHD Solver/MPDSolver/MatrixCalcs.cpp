void MulMatrix(double **MatA[5][5], double **MatB[5][5],
double **MatC[5][5], int Zone)
{
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            for(row=0; row<5; row++)
            {
                for(col=0; col<5; col++)
                {
                    MatC[row][col][j][k] = 0.0;
                    for(l=0; l<5; l++)
                    {
                        MatC[row][col][j][k] += MatA[row][l][j][k]*
                        MatB[l][col][j][k];
                    }
                }
            }
        }
    }
}

//-------------------------------------
void AddMatrix(double **MatA[5][5], double **MatB[5][5],
double **MatC[5][5], int Zone)
    {
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            for(row=0; row<5; row++)
            {
                for(col=0; col<5; col++)
                {
                    MatC[row][col][j][k] = MatA[row][col][j][k]
                    +MatB[row][col][j][k];
                }
            }
        }
    }
}

//-------------------------------------
void SubMatrix(double **MatA[5][5], double **MatB[5][5],
double **MatC[5][5], int Zone)
{
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            for(row=0; row<5; row++)
            {
                for(col=0; col<5; col++)
                {
                    MatC[row][col][j][k] = MatA[row][col][j][k]
                -MatB[row][col][j][k];
                }
            }
        }
    }
}

//-------------------------------------
void MatTimesVec(double **MatA[5][5], double **VecB[5],
double **VecC[5], int Zone)
{
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            for(row=0; row<5; row++)
            {
                VecC[row][j][k] = 0.0;
                for(col=0; col<5; col++)
                {
                    VecC[row][j][k] += MatA[row][col][j][k]
                    *VecB[col][j][k];
                }
            }
        }
    }
}

//-------------------------------------
void AddVector(double **VecA[5], double **VecB[5],
double **VecC[5], int Zone)
{
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            for(row=0; row<5; row++)
            {
                VecC[row][j][k] = VecA[row][j][k]+VecB[row][j][k];
            }
        }
    }
}

//-------------------------------------
void SubVector(double **VecA[5], double **VecB[5],
double **VecC[5], int Zone)
{
    for (j=Rmin[Zone]+1; j<=Rmax[Zone]; j++)
    {
        for (k=Zmin[Zone]+1; k<=Zmax[Zone]; k++)
        {
            for(row=0; row<5; row++)
            {
                VecC[row][j][k] = VecA[row][j][k]-VecB[row][j][k];
            }
        }
    }
}
//-------------------------------------
