// A Dynamic Programming based C++ program to find minimum
// number operations to convert str1 to str2
#include "mex.h"
#include <stack>

using namespace std;
// Utility function to find minimum of three numbers
int min(int x, int y, int z)
{
    if (x < y && x < z) return x;
    else if (y < z) return y;
    else return z;
}
void print2d(double *dp, size_t m, size_t n)
{   
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; j++)
        {
            mexPrintf("%2.0f ",dp[i+j*m]);
        }
        if (i > 0) mexPrintf("\b\n");
        else mexPrintf("\n");
    }
    mexPrintf("\n");
}
void print2d_int(int *dp, size_t m, size_t n)
{   
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; j++)
        {
            mexPrintf("%2d ",dp[i+j*m]);
        }
        if (i > 0) mexPrintf("\b\n");
        else mexPrintf("\n");
    }
    mexPrintf("\n");
}
void show_stack(stack <int> st)
{
    stack <int> tmp_st;
    while (st.size()>0)
    {
        tmp_st.push(st.top());
        mexPrintf("%4d\t",st.top());
        st.pop();
    }
    mexPrintf("\n");
    while (tmp_st.size()>0)
    {
        st.push(tmp_st.top());
        tmp_st.pop();
    }
}
stack <int> backtrace(size_t m, size_t n, double* dp, char *str1, char *str2)
{
    // Create a table to store results of subproblems
    stack <int> st,opts;
    int* visited = new int[m*n]();
    st.push(m*n-1);
    visited[m*n-1] = 1;
    while (st.size()>0)
    {
        //print2d_int(visited,m,n);
        //mexPrintf("Operations:\t");
        //show_stack(opts);
        int pos = st.top();
        st.pop();
        int i = pos % m;
        int j = (int)(pos / m);
        //mexPrintf("i = %d\tj = %d\n",i,j);
        if (i == 0 && j == 0) continue;
        //Row one Insertion
        if (i == 0) 
        {
            if (visited[i + m*(j-1)] == 0)
            {
                st.push(i + m*(j-1));
                visited[i + m*(j-1)] = 1;
            }
            opts.push(i + m*j + 3*m*n);
            continue;
        }
        //Col one Deletion
        if (j == 0) 
        {
            if (visited[i-1 + m*j] == 0)
            {
                st.push(i-1 + m*j);
                visited[i-1 + m*j] = 1;
            }
            opts.push(i + m*j + 2*m*n);
            continue;
        }
        //No change
        if (dp[i + m*j] == dp[i-1 + m*(j-1)] && visited[i-1 + m*(j-1)] == 0 && str1[i-1] == str2[j-1])
        {
            st.push(i-1 + m*(j-1));
            visited[i-1 + m*(j-1)] = 1;
            continue;
        }
        //Mutation
        if (dp[i + m*j] == (1 + dp[i-1 + m*(j-1)]))
        {
            if (visited[i-1 + m*(j-1)] == 0)
            {
                st.push(i-1 + m*(j-1));
                visited[i-1 + m*(j-1)] = 1;
            }
            opts.push(i + m*j + m*n);
        }
        //Deletion
        if (dp[i + m*j] == (1 + dp[i-1 + m*j]))
        {
            if (visited[i-1 + m*j] == 0)
            {
                st.push(i-1 + m*j);
                visited[i-1 + m*j] = 1;
            }
            opts.push(i + m*j + 2*m*n);
        }
        //Insertion
        if (dp[i + m*j] == (1 + dp[i + m*(j-1)]))
        {
            if (visited[i + m*(j-1)] == 0)
            {
                st.push(i + m*(j-1));
                visited[i + m*(j-1)] = 1;
            }
            opts.push(i + m*j + 3*m*n);
        }        
    }
    //show_stack(opts);
    delete visited;
    return opts;
}

void mexFunction(int nlhs, mxArray *plhs[], 
        int nrhs, const mxArray *prhs[])
{
    size_t m, n;
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    char *str1 = mxArrayToString(prhs[1]);
    char *str2 = mxArrayToString(prhs[2]);
    double* dp = mxGetPr(prhs[0]);
    stack <int> opts = backtrace(m,n,dp,str1,str2);
    plhs[0] = mxCreateNumericMatrix(opts.size(),1,mxINT32_CLASS, mxREAL);
    int* operation = (int *) mxGetData(plhs[0]);
    while (opts.size()>0)
    {
        operation[opts.size()-1] = opts.top();
        opts.pop();
    }
    return;
}
    