package bsim.physical.contact;

import javax.vecmath.Vector2d;
import javax.vecmath.Vector3d;

import bsim.BSim;
import bsim.capsule.BSimCapsuleBacterium;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import bsim.winter2021.Bacterium;
import bsim.winter2021.RectangleIntersection;

public class PhysicalContactBacterium extends Bacterium{

	public PhysicalContactBacterium(BSim sim, Vector3d px1, Vector3d px2) {
		super(sim, px1, px2);
		// TODO Auto-generated constructor stub
	}

	public PhysicalContactBacterium(BSim sim, Vector3d px1, Vector3d px2, long origin_id, long parent_id) {
		super(sim, px1, px2, origin_id, parent_id);
		// TODO Auto-generated constructor stub
	}

	private static double perimeter =0;
    // width of the bacterium (called radius because it also represents
    // how large the orbs are at each endpoint), all equal right now.
    public static double radius = 1.0/2.0;

//This is just a preliminary function, so it needs refinement
    public Vector3d NeighbourDistance(BSimCapsuleBacterium bacteria,BSimCapsuleBacterium neighbour_bac) {

        Vector3d dist = new Vector3d();
        dist.sub(bacteria.position, neighbour_bac.position); // position = midpoint of cell
      
        return dist;
    }
	

    public boolean isColliding(BSimCapsuleBacterium neighbour_bac, double r) {
        double bac_length = Math.sqrt(Math.pow(this.x2.x - this.x1.x, 2) + Math.pow(this.x2.y - this.x1.y, 2));
        double neighbour_bac_length = Math.sqrt(Math.pow(neighbour_bac.x2.x - neighbour_bac.x1.x, 2)
                + Math.pow(neighbour_bac.x2.y - neighbour_bac.x1.y, 2));
        double bac_angle = Math.atan((this.x2.y - this.x1.y) / (this.x2.x - this.x1.x));
        double neighbour_bac_angle = Math.atan((neighbour_bac.x2.y - neighbour_bac.x1.y) / (neighbour_bac.x2.x - neighbour_bac.x1.x));

        Vector2d[] rectangle_bac = RectangleIntersection.rectangle_vertices(this.position.x, this.position.y,
                bac_length + 2 * (radius + r), 2 * (radius + r), bac_angle);
        
        Vector2d[] rectangle_neighbour = RectangleIntersection.rectangle_vertices(neighbour_bac.position.x,
                neighbour_bac.position.y, neighbour_bac_length + 2 * (neighbour_bac.radius + r),
                2 * (neighbour_bac.radius + r), neighbour_bac_angle);
        
        
        if (RectangleIntersection.rectangle_intersection_perimeter(rectangle_bac, rectangle_neighbour)>perimeter) {      	
        	return true;
		}else {
			return false;
		}
    }    
    
    
    public boolean isCollidingFinal(BSimCapsuleBacterium neighbour_bac, double r,double area) {
        double bac_length = Math.sqrt(Math.pow(this.x2.x - this.x1.x, 2) + Math.pow(this.x2.y - this.x1.y, 2));
        double neighbour_bac_length = Math.sqrt(Math.pow(neighbour_bac.x2.x - neighbour_bac.x1.x, 2)
                + Math.pow(neighbour_bac.x2.y - neighbour_bac.x1.y, 2));
        double bac_angle = Math.atan((this.x2.y - this.x1.y) / (this.x2.x - this.x1.x));
        double neighbour_bac_angle = Math.atan((neighbour_bac.x2.y - neighbour_bac.x1.y) / (neighbour_bac.x2.x - neighbour_bac.x1.x));

        Vector2d[] rectangle_bac = RectangleIntersection.rectangle_vertices(this.position.x, this.position.y,
                bac_length + 2 * (radius + r), 2 * (radius + r), bac_angle);
        
        Vector2d[] rectangle_neighbour = RectangleIntersection.rectangle_vertices(neighbour_bac.position.x,
                neighbour_bac.position.y, neighbour_bac_length + 2 * (neighbour_bac.radius + r),
                2 * (neighbour_bac.radius + r), neighbour_bac_angle);
        
        
        if (rectangle_intersection_SutherlandHodgman(rectangle_bac, rectangle_neighbour)>area) {      	        	
        	return true;
		}else {
			return false;
		}
    }      
    
    
    
    private double rectangle_intersection_SutherlandHodgman(Vector2d[] r_1, Vector2d[] r_2) {
    	//Sutherland-Hodgman polygon clipping: https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
        List<double[]> subject, clipper, intersection;
    	List<double[]> subjPoints= new ArrayList<>();
        List<double[]> clipPoints = new ArrayList<>();
        
        //convert vector2d to list
        //converting counterclockwise vertices to clockwise vertices.
        for (int i = r_2.length-1; i >=0 ; i--) {
        	subjPoints.add(new double[]{r_2[i].x, r_2[i].y});
        }
        for (int i = r_1.length-1; i >=0 ; i--) {
        	clipPoints.add(new double[]{r_1[i].x, r_1[i].y});
        }
        
        subject = subjPoints;
        intersection  = subject; 
        clipper = clipPoints;          
        int len = clipper.size();        
 
        for (int i = 0; i < len; i++) {
 
            int len2 = intersection.size();
            List<double[]> input = intersection;
            intersection = new ArrayList<>(len2);
 
            double[] prev_clippper_vertex = clipper.get((i + len - 1) % len);
            double[] current_clippper_vertex = clipper.get(i);
 
            for (int j = 0; j < len2; j++) {
 
                double[] prev_input_vertex = input.get((j + len2 - 1) % len2);
                double[] current_input_vertex = input.get(j);
 
                if (isInside(prev_clippper_vertex, current_clippper_vertex, current_input_vertex)) {
                    if (!isInside(prev_clippper_vertex, current_clippper_vertex, prev_input_vertex))
                        intersection.add(intersectionPoint(prev_clippper_vertex, current_clippper_vertex, prev_input_vertex, current_input_vertex));
                    intersection.add(current_input_vertex);
                } else if (isInside(prev_clippper_vertex, current_clippper_vertex, prev_input_vertex))
                    intersection.add(intersectionPoint(prev_clippper_vertex, current_clippper_vertex, prev_input_vertex, current_input_vertex));
            }
        }
        return intersectionArea(intersection);
    }
    
    private boolean isInside(double[] a, double[] b, double[] c) {
   	 return (b[0]-a[0]) * (c[1]-a[1]) < (b[1]-a[1]) * (c[0]-a[0]);
   }
 
    private double[] intersectionPoint(double[] a, double[] b, double[] p, double[] q) {
        double num_x = (a[0]*b[1] - a[1]*b[0]) * (p[0]-q[0]) -
                (a[0]-b[0]) * (p[0]*q[1] - p[1]*q[0]);
      double den_x = (a[0]-b[0]) * (p[1]-q[1]) - (a[1]-b[1]) * (p[0]-q[0]);

      double num_y = (a[0]*b[1] - a[1]*b[0]) * (p[1]-q[1]) -
              (a[1]-b[1]) * (p[0]*q[1] - p[1]*q[0]);
      double den_y = (a[0]-b[0]) * (p[1]-q[1]) - (a[1]-b[1]) * (p[0]-q[0]);
      

      double x=num_x/den_x;
      double y=num_y/den_y;
      
      return new double[]{x, y};
    }
    
    private double intersectionArea(List<double[]> intersection) {
    	//I used Shoelace formula: https://en.wikipedia.org/wiki/Shoelace_formula
    	
    	double area=0;
    	
    	for (int i = 0; i < intersection.size()-1; i++) {
    		//x1.y2+x2.y3+x3.y4+....
    		area=area+(intersection.get(i)[0]*intersection.get(i+1)[1]);
    		//x2.y1-x3.y2-x4.y3-....
    		area=area-(intersection.get(i+1)[0]*intersection.get(i)[1]);
    		
    		if(i==intersection.size()-2) {
        		//xn.y1
        		area=area+(intersection.get(intersection.size()-1)[0]*intersection.get(0)[1]);
        		//x1.yn
        		area=area-(intersection.get(0)[0]*intersection.get(intersection.size()-1)[1]);    			
    		}
    	}
		return 0.5*Math.abs(area);
    }    
    

}
