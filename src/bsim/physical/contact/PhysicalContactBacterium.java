package bsim.physical.contact;

import javax.vecmath.Vector2d;
import javax.vecmath.Vector3d;

import bsim.capsule.BSimCapsuleBacterium;

import java.util.ArrayList;
import java.util.Arrays;

import bsim.winter2021.RectangleIntersection;

public class PhysicalContactBacterium {

    // coordinate vectors for the endpoints of the bacterium, initilaized at 0,0,0
    public Vector3d x1 = new Vector3d(0,0,0);
    public Vector3d x2 = new Vector3d(0,0,0);

    // center of the bacterium. note that position = (x2-x1)/2, initialize
    public Vector3d position = new Vector3d(0,0,0);
    
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
	

    public static boolean isColliding(BSimCapsuleBacterium bacteria,BSimCapsuleBacterium neighbour_bac, double r) {
        double bac_length = Math.sqrt(Math.pow(bacteria.x2.x - bacteria.x1.x, 2) + Math.pow(bacteria.x2.y - bacteria.x1.y, 2));
        double neighbour_bac_length = Math.sqrt(Math.pow(neighbour_bac.x2.x - neighbour_bac.x1.x, 2)
                + Math.pow(neighbour_bac.x2.y - neighbour_bac.x1.y, 2));
        double bac_angle = Math.atan((bacteria.x2.y - bacteria.x1.y) / (bacteria.x2.x - bacteria.x1.x));
        double neighbour_bac_angle = Math.atan((neighbour_bac.x2.y - neighbour_bac.x1.y) / (neighbour_bac.x2.x - neighbour_bac.x1.x));

        Vector2d[] rectangle_bac = RectangleIntersection.rectangle_vertices(bacteria.position.x, bacteria.position.y,
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

}
