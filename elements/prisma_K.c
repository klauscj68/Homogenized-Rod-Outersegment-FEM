 	#include "mex.h"
	#include "math.h"

	void prisma_K(double X[], double Y[], double Z[], double K[][6])
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
	  
	  
      t1 = Z[0]-Z[1];
      t2 = X[2]*X[2];
      t3 = Y[1]*Y[1];
      t4 = Y[2]*Y[2];
      t5 = X[1]*X[1];
      t6 = Y[1]*Y[2];
      t8 = X[1]*X[2];
      t19 = 1/(-X[0]*Y[2]-X[1]*Y[0]+X[1]*Y[2]+X[0]*Y[1]+X[2]*Y[0]-X[2]*Y[1]);
      t20 = t1*(t2+t3+t4+t5-2.0*t6-2.0*t8)*t19;
      t21 = t20/6.0;
      t22 = Y[0]*Y[2];
      t23 = Y[0]*Y[1];
      t24 = X[0]*X[2];
      t25 = X[0]*X[1];
      t28 = t1*(t4+t2-t22+t23-t6-t24+t25-t8)*t19;
      t29 = t28/6.0;
      t32 = t1*(-t5-t3-t22+t6+t23-t24+t8+t25)*t19;
      t33 = t32/6.0;
      t34 = t20/12.0;
      t35 = t28/12.0;
      t36 = t32/12.0;
      t37 = X[0]*X[0];
      t38 = Y[0]*Y[0];
      t43 = t1*(t37+t4+t38+t2-2.0*t22-2.0*t24)*t19;
      t44 = t43/6.0;
      t47 = t1*(t37+t38-t22+t6-t23-t24+t8-t25)*t19;
      t48 = t47/6.0;
      t49 = t43/12.0;
      t50 = t47/12.0;
      t55 = t1*(t37+t38+t5+t3-2.0*t23-2.0*t25)*t19;
      t56 = t55/6.0;
      t57 = t55/12.0;
      K[0][0] = -t21;
      K[0][1] = t29;
      K[0][2] = -t33;
      K[0][3] = -t34;
      K[0][4] = t35;
      K[0][5] = -t36;
      K[1][0] = t29;
      K[1][1] = -t44;
      K[1][2] = t48;
      K[1][3] = t35;
      K[1][4] = -t49;
      K[1][5] = t50;
      K[2][0] = -t33;
      K[2][1] = t48;
      K[2][2] = -t56;
      K[2][3] = -t36;
      K[2][4] = t50;
      K[2][5] = -t57;
      K[3][0] = -t34;
      K[3][1] = t35;
      K[3][2] = -t36;
      K[3][3] = -t21;
      K[3][4] = t29;
      K[3][5] = -t33;
      K[4][0] = t35;
      K[4][1] = -t49;
      K[4][2] = t50;
      K[4][3] = t29;
      K[4][4] = -t44;
      K[4][5] = t48;
      K[5][0] = -t36;
      K[5][1] = t50;
      K[5][2] = -t57;
      K[5][3] = -t33;
      K[5][4] = t48;
      K[5][5] = -t56;
    }


	/* The gateway routine */
    void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

    {
        double *X, *Y, *Z, *K;
   
		/* Check for proper number of arguments */

		if (nrhs != 3) {
			mexErrMsgTxt("Three input arguments required.");
		}
		else if (nlhs != 1) {
			mexErrMsgTxt("One output argument required.");
		}

		X = mxGetPr(prhs[0]);
        Y = mxGetPr(prhs[1]);
		Z = mxGetPr(prhs[2]);

		
		/* Create a matrix for the return argument */
		plhs[0] = mxCreateDoubleMatrix(6,6,mxREAL);

		/* Assign pointers to each input and output. */
		K = mxGetPr(plhs[0]);

		/* Call the C subroutine */
        prisma_K(X,Y,Z,K);
	}
