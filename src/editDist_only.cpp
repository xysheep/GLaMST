// A Dynamic Programming based C++ program to find minimum
// number operations to convert str1 to str2
#include "mex.h"

// Utility function to find minimum of three numbers
int min(int x, int y, int z)
{
    if (x < y && x < z) return x;
    else if (y < z) return y;
    else return z;
}
void print2d(double *dp, size_t m, size_t n)
{   
    for (int i=0; i<=m; i++)
    {
        for (int j=0; j<=n; j++)
        {
            mexPrintf("%2.0f ",dp[i+j*(m+1)]);
        }
        if (i > 0) mexPrintf("\b\n");
        else mexPrintf("\n");
    }
    mexPrintf("\n");
}
int editDistDP(char *str1, char *str2, size_t m, size_t n, double* dp)
{
    // Create a table to store results of subproblems
    
    
    for (int i = 0;i<= (m+1)*(n+1)-1;i++)dp[i]= 9999;
//     print2d(dp, m, n);
    // Fill d[][] in bottom up manner
    for (int i=0; i<=m; i++)
    {
        for (int j=0; j<=n; j++)
        {
            // If first string is empty, only option is to
            // isnert all characters of second string
            if (i==0)
                dp[i+j*(m+1)] = j;  // Min. operations = j
 
            // If second string is empty, only option is to
            // remove all characters of second string
            else if (j==0)
                dp[i+j*(m+1)] = i; // Min. operations = i
 
            // If last characters are same, ignore last char
            // and recur for remaining string
            else if (str1[i-1] == str2[j-1])
                dp[i+j*(m+1)] = dp[(i-1)+(j-1)*(m+1)];
 
            // If last character are different, consider all
            // possibilities and find minimum
            else
                dp[i+j*(m+1)] = 1 + min(dp[i+(j-1)*(m+1)],  // Insert
                                   dp[(i-1)+j*(m+1)],  // Remove
                                   dp[(i-1)+(j-1)*(m+1)]); // Replace
//              print2d(dp, m, n);
        }
    }
    int dist = dp[(m+1)*(n+1)-1];
    return dist;
}

void mexFunction(int nlhs, mxArray *plhs[], 
        int nrhs, const mxArray *prhs[])
{
    size_t n1, n2;
    n1 = mxGetN(prhs[0]);
    n2 = mxGetN(prhs[1]);
    plhs[1] = mxCreateNumericMatrix(n1+1,n2+1,mxDOUBLE_CLASS, mxREAL);
    double* dp = (double *) mxGetData(plhs[1]);
    char *str1 = mxArrayToString(prhs[0]);
    char *str2 = mxArrayToString(prhs[1]);
    plhs[0] = mxCreateDoubleScalar(editDistDP(str1,str2,n1,n2,dp));
    mxFree(str1);
    mxFree(str2);
    return;
}
    