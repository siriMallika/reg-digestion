
import java.util.ArrayList;
import java.util.List;

import org.opencv.core.Core;
import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.Size;
import org.opencv.imgcodecs.Imgcodecs;
import org.opencv.imgproc.Imgproc;

class ImgReg{
	  public  void run() {
		// TODO Auto-generated method stub
			float resizeImg=0.5f;
		    //int num=256;
		    //double h1,h2;
		    int [] k= new int[31];
		   
		    //Mat t = Imgcodecs.imread(getClass().getResource("/eS001.0-31_Frame1.jpg").getPath());
		    //Mat s = Imgcodecs.imread(getClass().getResource("/Dynamic32-37min.0146.jpg").getPath());
		    //Mat s = Imgcodecs.imread(getClass().getResource("/static60min.jpg").getPath());
		    //Mat s = Imgcodecs.imread(getClass().getResource("/eS001.0-31_Frame31.jpg").getPath());
			Mat s = Imgcodecs.imread(getClass().getResource("/circle2.jpg").getPath());
		    Mat t = Imgcodecs.imread(getClass().getResource("/circle1.jpg").getPath());
		    Mat gS = new Mat(s.height(),s.width(), CvType.CV_8U);
		    Mat gT = new Mat(t.height(),t.width(), CvType.CV_8U);
		    Mat gSm = new Mat(s.height(),s.width(), CvType.CV_8U);
		    Mat gTm = new Mat(t.height(),t.width(), CvType.CV_8U);
		    Imgproc.cvtColor(s, gS, Imgproc.COLOR_RGB2GRAY);
		    Imgproc.cvtColor(t, gT, Imgproc.COLOR_RGB2GRAY);
		    
		    Size dsize0=gS.size();
		    dsize0.height=dsize0.height*resizeImg;
		    dsize0.width=dsize0.width*resizeImg;
		    Imgproc.resize(gS, gSm, dsize0);
		    Imgproc.resize(gT, gTm, dsize0);
		    List<Mat> Sa = new ArrayList<Mat>(3);
		    Core.split(gSm, Sa);
		    System.out.println("S Size= "+Sa.size());
		    Mat channelS = Sa.get(0);
		    System.out.println("S= "+channelS.get(0,1)[0]);
		    List<Mat> Ta = new ArrayList<Mat>(3);
		    Core.split(gTm, Ta);
		    System.out.println("T size= "+Ta.size());
		    Mat channelT = Ta.get(0);
		    System.out.println("T= "+channelT.get(0,1)[0]);
		   
		    
		    clMLrigidreg reg= new clMLrigidreg();
		    double maxt=16.0,maxr=16.0;
		    reg.init (gSm,gTm);
		    double r=reg.run_ncc(maxt, maxr, 1);
		    System.out.println("cost= "+r);
		    System.out.println("m_nP0 R= "+reg.m_nP[0]);
		    System.out.println("m_nP1 Tx= "+reg.m_nP[1]);
		    System.out.println("m_nP2 Ty= "+reg.m_nP[2]);
		    
		    int X=channelS.width();
			int Y=channelS.height();
			Mat img2 = new Mat(Y,X, CvType.CV_8U);
		       //reg.m_m=new double [Y][X][1];                 // target image
		    
			  for(int i=0;i<Y;i++){
				  for(int j=0;j<X;j++){
					 // System.out.print(i +" " +j);
					 double data=reg.m_m[i][j][0];
						img2.put(i, j, data);
					 
				  }
			  }
		    
		    String filename1 = "imgResult1-1.png";
		    Imgcodecs.imwrite(filename1, img2);
		    
		     img2 = new Mat(Y,X, CvType.CV_8U);
		       //reg.m_m=new double [Y][X][1];                 // target image
		    
			  for(int i=0;i<Y;i++){
				  for(int j=0;j<X;j++){
					 // System.out.print(i +" " +j);
					 double data=reg.m_t[i][j][0];
						img2.put(i, j, data);
					 
				  }
			  }
		    
		     filename1 = "imgTarget1-2.png";
		    Imgcodecs.imwrite(filename1, img2);
		    
		     img2 = new Mat(Y,X, CvType.CV_8U);
		       //reg.m_m=new double [Y][X][1];                 // target image
		    
			  for(int i=0;i<Y;i++){
				  for(int j=0;j<X;j++){
					 // System.out.print(i +" " +j);
					 double data=reg.m_s[i][j][0];
						img2.put(i, j, data);
					 
				  }
			  }
		    
		     filename1 = "imgSource1-2.png";
		    Imgcodecs.imwrite(filename1, img2);
		    
		    double[][][] movei=new double [Y][X][1];
		    movei=reg.compute_moving_image1(reg.m_m,reg.m_nP);
		    
		    img2 = new Mat(Y,X, CvType.CV_8U);
		       //reg.m_m=new double [Y][X][1];                 // target image
		    
			  for(int i=0;i<Y;i++){
				  for(int j=0;j<X;j++){
					 // System.out.print(i +" " +j);
					 double data=movei[i][j][0];
						img2.put(i, j, data);
					 
				  }
			  }
			  
			  filename1 = "imgResult2-2.png";
			    Imgcodecs.imwrite(filename1, img2);
//		    Mat yuvimg0 = new Mat(image0.height(),image0.width(), CvType.CV_8U);
//		    Mat smallPic0 = new Mat(image0.height(),image0.width(), CvType.CV_8U);
//		    
//		    Size dsize0=image0.size();
//		    dsize0.height=dsize0.height*resizeImg;
//		    dsize0.width=dsize0.width*resizeImg;
//		    Imgproc.resize(image0, smallPic0, dsize0);
//		    Imgproc.cvtColor(smallPic0, yuvimg0, Imgproc.COLOR_RGB2YUV);
//		    List<Mat> Yimg0 = new ArrayList<Mat>(3);
//		    Core.split(yuvimg0, Yimg0);
//
//		    Mat channel0 = Yimg0.get(0);
//		    Mat imgAvg = new Mat(channel0.height(),channel0.width(), CvType.CV_8U);
//		    Mat imgSD2 = new Mat(channel0.height(),channel0.width(), CvType.CV_8U);
//		    double[] data=new double[31];
//		    double[] a=new double[31];
//		    
//		    int X=channel0.width();
//			  int Y=channel0.height();
//			  
//
//	    String filename1 = "imgSD21.png";
//	    Imgcodecs.imwrite(filename1, imgSD2);
		   
		}
}

public class ImageReg {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		System.out.println("Hello, OpenCV");
	    // Load the native library.
	    System.loadLibrary(Core.NATIVE_LIBRARY_NAME);
	    new ImgReg().run();

	}

}
