 	#include "mex.h"
	#include "math.h"

	void rett_K(double S[], double Z[], double K[][4])
	{
	
    double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    double t10, t11, t12, t13, t14, t15, t16, t17, t18, t19;
    double t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
    double t30, t31, t32, t33, t34, t35, t36, t37, t38, t39;
    double t40, t41, t42, t43, t44, t45, t46, t47, t48, t49;
    double t50, t51, t52, t53, t54, t55, t56, t57, t58, t59;
    double t60, t61, t62, t63, t64, t65, t66, t67, t68, t69;
    double t70, t71, t72, t73, t74, t75, t76, t77, t78, t79;
    double t80, t81, t82, t83, t84, t85, t86, t87, t88, t89;
    double t90, t91, t92, t93, t94, t95, t96, t97, t98, t99;
	  
	  
      t1 = S[0]*S[0];
      t2 = t1*S[0];
      t3 = 3.0*t2;
      t4 = S[1]*S[1];
      t5 = S[1]*t4;
      t6 = Z[1]*Z[1];
      t7 = t6*S[0];
      t8 = 2.0*t7;
      t9 = Z[0]*Z[0];
      t10 = t9*S[1];
      t11 = 2.0*t10;
      t12 = t1*S[1];
      t13 = 5.0*t12;
      t14 = t6*S[1];
      t15 = 2.0*t14;
      t16 = t9*S[0];
      t17 = 2.0*t16;
      t18 = S[0]*t4;
      t19 = Z[0]*Z[1];
      t20 = t19*S[0];
      t21 = 4.0*t20;
      t22 = t19*S[1];
      t23 = 4.0*t22;
      t26 = 1/(S[0]-S[1]);
      t29 = 1/(Z[0]-Z[1]);
      t31 = (t3+t5+t8+t11-t13+t15+t17+t18-t21-t23)*t26*t29/12.0;
      t32 = S[1]+S[0];
      t34 = 2.0*S[0]*S[1];
      t40 = t26*t29;
      t42 = t32*(t1-t34+t4+4.0*t19-2.0*t6-2.0*t9)*t40/12.0;
      t43 = 2.0*t20;
      t44 = 2.0*t22;
      t48 = (t5+t3-t13+t18-t16-t10-t7-t14+t43+t44)*t26*t29/12.0;
      t53 = t32*(t1-t34+t6+t9+t4-2.0*t19)*t40/12.0;
      t54 = 3.0*t5;
      t55 = 5.0*t18;
      t59 = (t54+t2-t55+t12+t17+t11+t8+t15-t23-t21)*t26*t29/12.0;
      t63 = (t54+t2+t12-t55-t16-t10-t7-t14+t43+t44)*t26*t29/12.0;
      K[0][0] = t31;
      K[0][1] = t42;
      K[0][2] = -t48;
      K[0][3] = -t53;
      K[1][0] = t42;
      K[1][1] = t59;
      K[1][2] = -t53;
      K[1][3] = -t63;
      K[2][0] = -t48;
      K[2][1] = -t53;
      K[2][2] = t31;
      K[2][3] = t42;
      K[3][0] = -t53;
      K[3][1] = -t63;
      K[3][2] = t42;
      K[3][3] = t59;
    }


	/* The gateway routine */
    void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

    {
        double *S, *Z, *K;
   
		/* Check for proper number of arguments */

		if (nrhs != 2) {
			mexErrMsgTxt("Two input arguments required.");
		}
		else if (nlhs != 1) {
			mexErrMsgTxt("One output argument required.");
		}

		S = mxGetPr(prhs[0]);
        Z = mxGetPr(prhs[1]);

		
		/* Create a matrix for the return argument */
		plhs[0] = mxCreateDoubleMatrix(4,4,mxREAL);

		/* Assign pointers to each input and output. */
		K = mxGetPr(plhs[0]);

		/* Call the C subroutine */
	    rett_K(S, Z, K);
	}
