package bsim.physical.contact;

import javax.vecmath.Vector2d;
import javax.vecmath.Vector3d;

import bsim.BSim;
import bsim.capsule.BSimCapsuleBacterium;

import java.util.ArrayList;
import java.util.Arrays;

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

}
