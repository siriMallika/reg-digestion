import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.opencv.core.Core;
import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.Size;
import org.opencv.imgcodecs.Imgcodecs;
import org.opencv.imgproc.Imgproc;


public class clMLrigidreg {
	double [][][]    m_s;                // source image
    double [][][]   m_t;                // target image
    double [][][]   m_m;                // moving image
    double          m_rt, m_rr;         // translation and rotation factor
    int             m_cx, m_cy, m_nc;   // image details
    double [][]    m_RT = new double [3][3];   // rigid body transformation matrix
    double [][]    m_iT = new double [3][3];   // inversed rigid body transformation matrix
    double []       m_nP = new double [3];      // normalized optimizing parameters [R, Tx, Ty]
    // results
  //  double          RotationRad { get { return m_nP [0]*m_rr; } }
  //  double          RotationDeg { get { return m_nP [0]*m_rr*180.0/Math.PI; } }
  // double []       Translation { get { return new double [] { m_nP [1]*m_rt, m_nP [2]*m_rt }; } }

  //  public void init (Image <Gray, double> s, Image <Gray, double> t)
    public void init (Mat s,Mat t)
    {
        if (s == null || t == null) return;

        //plus+++

	    List<Mat> Sa = new ArrayList<Mat>(3);
	    Core.split(s, Sa);
	    System.out.println("S Size= "+Sa.size());
	    Mat channelS = Sa.get(0);
	    System.out.println("S= "+channelS.get(0,0)[0]);
	    
	    List<Mat> Ta = new ArrayList<Mat>(3);
	    Core.split(t, Ta);
	    System.out.println("T size= "+Ta.size());
	    Mat channelT = Ta.get(0);
	    System.out.println("T= "+channelT.get(0,1)[0]);
	    
	    int X=channelS.width();
		int Y=channelS.height();
	
	      m_s=new double [Y][X][1];                // source image
	       m_t=new double [Y][X][1];                 // target image
	    
		  for(int i=0;i<Y;i++){
			  for(int j=0;j<X;j++){
				 // System.out.print(i +" " +j);
				  m_s[i][j][0] = channelS.get(i,j)[0];
				  m_t[i][j][0] = channelT.get(i,j)[0];
			  }
		  }

		  
        // initialize images (source s, target t and moving m)
	   // m_s=Sa.get(0);
	   // m_t=Ta.toArray(m_t);
        //m_s     = (double [][][]) Sa.toArray(m_s);
        //m_s     = (double [][][]) s.Data.Clone ();
        //m_t     = (double [][][]) t.Clone ();
        //m_t     = (double [][][]) t.Data.Clone ();
		  m_cy    = Y;    // height
	        m_cx    = X;    // width
	        m_nc    = 1;    // color components
//        m_cy    = m_s[0].length;    // height
//        m_cx    = m_s[1].length;    // width
//        m_nc    = m_s[2].length;    // color components
        m_m     = new double [m_cy][ m_cx][ m_nc];
        // initialize scaling factors
        m_rt    = 1.0;  // translation scale
        m_rr    = 1.0;  // rotation scale
        init_Tx ();
    }

    // return inverse of rigid body transformation T in homogeneous coordinates
    private double [][]  inverseTX (double [][] T)
    {
        double [] [] ret = new double [] []{ { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };   // identity matrix

        ret [0] [0] = T [0] [0];  ret [0] [1] = T [1][ 0];  // R = R^T
        ret [1][ 0] = T [0][ 1];  ret [1][ 1] = T [1][ 1];

        ret [0][ 2] = -(ret [0][ 0]*T [0][ 2] + ret [0][ 1]*T [1][ 2]);  // T = -(R^T)*T
        ret [1][ 2] = -(ret [1][ 0]*T [0][ 2] + ret [1][ 1]*T [1][ 2]);

        return ret;
    }

    // covert normalized parameter to 2D rigid body transformation 
    private double [][ ] p2T (double [] p)
    {
        double [][ ] ret = new double [][ ] { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };   // identity matrix
        double      tx, ty, cs, sn; // cs = cos (theta), sn = sin (theta); theta = p[0]*m_rr

        cs = Math.cos (p [0]*m_rr);
        sn = Math.sin (p [0]*m_rr);
        tx = p [1]*m_rt;
        ty = p [2]*m_rt;
        // construct rigid body tranx
        ret [0][ 0] =  cs;   ret [0][ 1] = -sn;   ret [0][ 2] = tx;
        ret [1][ 0] =  sn;   ret [1][ 1] =  cs;   ret [1][ 2] = ty;

        return ret;
    }

    // initialize transformation and related aux variables to indentity
    private void init_Tx ()
    {
        m_RT = new double [][ ] {{ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 }};   // identity matrix
        m_nP = new double [] { 0.0, 0.0, 0.0 };     // zero translation and rotation

        m_RT = p2T (m_nP);
        m_iT = inverseTX (m_RT);
    }

    
    
    // rigid body transformation of a pixel (sx, sy) using transformation T
    // the result is returned to (dx, dy)
   // private void Tx (double [][ ] T, double sx, double sy, out double dx, out double dy)
    private double Txdx (double [][ ] T, double sx, double sy)
    {
    	//plus+++
    	//Swapping p =new Swapping();
    	//
        return (sx*T [0][ 0] + sy*T [0][ 1] + T [0][ 2]);
       // p.dy = sx*T [1][ 0] + sy*T [1][ 1] + T [1][ 2];
    }
    private double Txdy (double [][ ] T, double sx, double sy)
    {
    	//plus+++
    	//Swapping p =new Swapping();
    	//
        //dx = sx*T [0][ 0] + sy*T [0][ 1] + T [0][ 2];
    	return ( sx*T [1][ 0] + sy*T [1][ 1] + T [1][ 2]);
    }

    // fetch image pixel (component c) at a given (double) coordinate (dx, dy)
    private double fetch_source (double dx, double dy, int c )
    {
    	c = 0;
        int     ix = (int) (dx + 0.5), iy = (int) (dy + 0.5);

        // boundary assertion
        if (ix < 0 || ix > m_cx - 1 || iy < 0 || iy > m_cy - 1)
            return 0;
        // neearest neighbor interpolation
        return m_s [iy][ ix][ c];
    }

    // transform the source image into a moving image based on parameter p
    // use inverseTX to find source coordinates in m_s from a given destination ones in mvi
    private double compute_moving_image (double [][][ ] mvi, double [] p)
    {
        double [][ ] iT = inverseTX (p2T (p));   // calculate inverse rigid body from parameter p
        int         x, y;
        double      dx, dy, avg = 0.0;
      //plus+++
    	//Swapping sp =new Swapping();
    	//

        for (y = 0; y < m_cy; y ++) {           // loop to each pixel in moving image
            for (x = 0; x < m_cx; x ++) {
                //Tx (iT, x, y, sp.dx, sp.dy);  // get the corresponding coordinate in the source image
            	dx=Txdx (iT, x, y);
            	dy=Txdy (iT, x, y);
                mvi [y][ x][ 0] = fetch_source (dx, dy, 0);
                avg = avg + mvi [y][ x][ 0];      // running average for validation
            }
        }
        return (avg/(m_cx*m_cy));
    }
    
 // transform the source image into a moving image based on parameter p
    // use inverseTX to find source coordinates in m_s from a given destination ones in mvi
    double[][][] compute_moving_image1 (double [][][ ] mvi, double [] p)
    {
        double [][ ] iT = inverseTX (p2T (p));   // calculate inverse rigid body from parameter p
        int         x, y;
        double      dx, dy, avg = 0.0;
      
        for (y = 0; y < m_cy; y ++) {           // loop to each pixel in moving image
            for (x = 0; x < m_cx; x ++) {
                //Tx (iT, x, y, sp.dx, sp.dy);  // get the corresponding coordinate in the source image
            	dx=Txdx (iT, x, y);
            	dy=Txdy (iT, x, y);
                mvi [y][ x][ 0] = fetch_source (dx, dy, 0);
                //avg = avg + mvi [y][ x][ 0];      // running average for validation
            }
        }
        return mvi;
    }

    // compute inverted normalized correlation coefficient between source (s) and target (t)
    public double ncc (double [][][ ] s, double [][][] t)
    {
        int     x, y, omega = m_cx*m_cy;
        //double  ncorr = float.MaxValue;
        double  ncorr =Float.MAX_VALUE;
        double  means = 0.0, meant = 0.0;
        double  ivals = 0.0, ivalt = 0.0;
        double  s1 = 0.0, s2 = 0.0, s3 = 0.0;

        // mean of both images
        for (y = 0; y < m_cy; y++)
            for (x = 0; x < m_cx; x++) {
                means = means + s [y][ x][ 0];
                meant = meant + t [y][ x][ 0];
            }
        means = means/ omega;
        meant = meant/ omega;

        // compute (auto, cross) correlation
        for (y = 0; y < m_cy; y++)
            for (x = 0; x < m_cx; x++) {
                ivals = (s [y][ x][ 0] - means);
                ivalt = (t [y][ x][ 0] - meant);

                s1 = s1 + ivals*ivalt;
                s2 = s2 + ivals*ivals;
                s3 = s3 + ivalt*ivalt;
            }

        ncorr = (s2*s3 > 0.0) ? s1/Math.sqrt (s2*s3) : 0.0;
        return 1.0 - Math.abs (ncorr);
    }

    // cost function with respect to parameter p
    public double cost (double [] p)
    {
        double  avg, sim;
        // convert parameter p to transformation and prepare the respective moving image
        avg = compute_moving_image (m_m, p);
        // return normalized cross-correlation (ncc)
        sim = ncc (m_m, m_t);

        return sim;
    }

    // random next plausible value of parameter vector (-1.0 to 1.0)
    public void rand_nextp (Random rnd, double [] n)
    {
        for (int i = 0; i < n.length; i ++) 
            n [i] = 2.0*rnd.nextDouble () - 1.0f;
    }
    
    // optimize translation and rotation parameter with maximum pixels and degrees, respectively
    // the default maximum iterations are set to 1024 
    public double run_ncc (double maxt, double maxr, int iterations)
    {
    	iterations = 1024;
        double []   curr    = new double [3];
        double []   next    = new double [3];

        //the probability
        double      proba, alpha = 1.0 - (1.0/iterations), temperature = 1.0, epsilon = 0.001, delta;
        double      ccost   = cost (curr);
        double      gcost   = ccost;
        Random      rnd     = new Random ();

        m_rt = maxt;
        m_rr = maxr*Math.PI/180.0;
        init_Tx ();

        while (temperature > epsilon) {
            for (int i = 0; i < iterations; i ++) {
                rand_nextp (rnd, curr);         // get the next random permutation of config
                delta = cost (curr) - ccost;    // compute the distance of the new permuted config

                if (delta < 0) {                // if the new distance is better accept it and assign it
                	System.arraycopy(curr, 0, next,0, curr.length);
                	//Array.Copy (curr, next, curr.length);
                    ccost   = delta + ccost;
                } else {
                    proba = rnd.nextDouble ();
                    if (proba < Math.exp (-delta / temperature)) {
                    	System.arraycopy(curr, 0, next,0, curr.length);
                        //Array.Copy (curr, next, curr.length);
                        ccost   = delta + ccost;
                    }
                }
                temperature *= alpha;
            }
            // update value at every temperature
            if (ccost < gcost) {
            	System.arraycopy(next, 0, m_nP,0, next.length);
                //Array.Copy (next, m_nP, next.length);
                gcost = ccost;
            }
        }

        return gcost;
    }
}

////plus++++
//class Swapping{
//	private double dx;
//	private double dy;
//	public double getDx() {
//		return dx;
//	}
//	public void setDx(double dx) {
//		this.dx = dx;
//	}
//	public double getDy() {
//		return dy;
//	}
//	public void setDy(double dy) {
//		this.dy = dy;
//	}
//	
//}
