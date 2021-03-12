// NAME: Zhiguo(Frank) Zhang
// ID: 260550226
package comp559.project;

import javax.vecmath.Vector3d;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;

/**
 * Spring class for 599 assignment 1
 * @author kry
 */
public class Spring {

    Particle p1 = null;
    Particle p2 = null;
    
    boolean boundary = false;
    /** Spring stiffness, sometimes written k_s in equations */
    public static double k = 1;
    /** Spring damping (along spring direction), sometimes written k_d in equations */
    public static double c = 1;
    /** Rest length of this spring */
    double l0 = 0;
    
    /**
     * Creates a spring between two particles
     * @param p1
     * @param p2
     */
    public Spring( Particle p1, Particle p2) {
        this.p1 = p1;
        this.p2 = p2;
        recomputeRestLength();
        p1.springs.add(this);
        p2.springs.add(this);
    }
    
    public Spring( Particle p1, Particle p2, boolean boundary) {
        this.p1 = p1;
        this.p2 = p2;
        this.boundary = boundary;
        recomputeRestLength();
        p1.springs.add(this);
        p2.springs.add(this);
    }
    
    /**
     * Computes and sets the rest length based on the original position of the two particles 
     */
    public void recomputeRestLength() {
        l0 = p1.p0.distance( p2.p0 );
    }
    
    /**
     * Applies the spring force by adding a force to each particle
     */
    public void apply() {
        // TODO: Objective 1, FINISH THIS CODE!
    		Vector3d v1 = new Vector3d();
    		v1.set(p1.v);
    		v1.sub(p2.v);
    		Vector3d l = new Vector3d();
    		l.set(p1.p);
    		l.sub(p2.p);
        double scale = - (k * (l.length() - l0) + c * v1.dot(l) / l.length()) / l.length();
        l.scale(scale);
        p1.addForce(l);
        l.negate();
        p2.addForce(l);
    }
       
    /**
     * Computes the force and adds it to the appropriate components of the force vector.
     * (This function is something you might use for a backward Euler integrator)
     * @param f
     */
    public void addForce(Vector f) {
        // TODO: Objective 8, FINISH THIS CODE for backward Euler method (probably very simlar to what you did above)
    		Vector3d v1 = new Vector3d();
		v1.set(p1.v);
		v1.sub(p2.v);
		Vector3d l = new Vector3d();
		l.set(p1.p);
		l.sub(p2.p);
	    double scale = - (k * (l.length() - l0) + c * v1.dot(l) / l.length()) / l.length();
	    l.scale(scale);
	    f.add(3*p1.index, l.x);
	    f.add(3*p1.index+1, l.y);
	    f.add(3*p1.index+2, l.z);
	    l.negate();
	    f.add(3*p2.index, l.x);
	    f.add(3*p2.index+1, l.y);
	    f.add(3*p2.index+2, l.z);
    }
    
    /**
     * Adds this springs contribution to the stiffness matrix
     * @param dfdx
     */
    public void addDfdx( Matrix dfdx ) {
        // TODO: Objective 8, FINISH THIS CODE... necessary for backward euler integration
	    DenseVector xhat = new DenseVector(3);
	    	DenseMatrix xxt = new DenseMatrix(3,3);
	    	DenseMatrix A = new DenseMatrix(3,3);
	    	
	    	// l = p1 - p2
	    	Vector3d l = new Vector3d();
		l.set(p1.p);
		l.sub(p2.p);
		double len = l.length();
		l.normalize();
		
		// xhat = lhat 
		xhat.set(0,l.x);
		xhat.set(1,l.y);
		xhat.set(2,l.z);
		
		// xxt = x * xt
	    	xxt.rank1(xhat);
	    	
	    	// dfdx := - k * (( 1 - l0/len ) * (I - x * xt) + x * xt)
	    	// A := ( 1 - l0/len ) * (I - x * xt)
	    	// xxt := - k * ( A + xxt )
	    	A.set(xxt);
	    	A.scale(-1);
	    	A.add(0, 0, 1.0);
	    	A.add(1, 1, 1.0);
	    	A.add(2, 2, 1.0);
	    	A.scale(1 - l0/len);
	    	xxt.add(A);
	    	xxt.scale(-k);
	    	
	    	int p1x = 3 * p1.index;
	    	int p1y = 3 * p1.index + 1;
	    	int p1z = 3 * p1.index + 2;
	    	int p2x = 3 * p2.index;
	    	int p2y = 3 * p2.index + 1;
	    	int p2z = 3 * p2.index + 2;
	    	
	    	dfdx.add(p1x, p1x, xxt.get(0, 0));
	    	dfdx.add(p1x, p1y, xxt.get(0, 1));
	    	dfdx.add(p1x, p1z, xxt.get(0, 2));
	    	dfdx.add(p1y, p1x, xxt.get(1, 0));
	    	dfdx.add(p1y, p1y, xxt.get(1, 1));
	    	dfdx.add(p1y, p1z, xxt.get(1, 2));
	    	dfdx.add(p1z, p1x, xxt.get(2, 0));
	    	dfdx.add(p1z, p1y, xxt.get(2, 1));
	    	dfdx.add(p1z, p1z, xxt.get(2, 2));
	    	
	    	dfdx.add(p2x, p2x, xxt.get(0, 0));
	    	dfdx.add(p2x, p2y, xxt.get(0, 1));
	    	dfdx.add(p2x, p2z, xxt.get(0, 2));
	    	dfdx.add(p2y, p2x, xxt.get(1, 0));
	    	dfdx.add(p2y, p2y, xxt.get(1, 1));
	    	dfdx.add(p2y, p2z, xxt.get(1, 2));
	    	dfdx.add(p2z, p2x, xxt.get(2, 0));
	    	dfdx.add(p2z, p2y, xxt.get(2, 1));
	    	dfdx.add(p2z, p2z, xxt.get(2, 2));

	    	// dfdx_i = - dfdx_j
	    	xxt.scale(-1);
	    	dfdx.add(p1x, p2x, xxt.get(0, 0));
	    	dfdx.add(p1x, p2y, xxt.get(0, 1));
	    	dfdx.add(p1x, p2z, xxt.get(0, 2));
	    	dfdx.add(p1y, p2x, xxt.get(1, 0));
	    	dfdx.add(p1y, p2y, xxt.get(1, 1));
	    	dfdx.add(p1y, p2z, xxt.get(1, 2));
	    	dfdx.add(p1z, p2x, xxt.get(2, 0));
	    	dfdx.add(p1z, p2y, xxt.get(2, 1));
	    	dfdx.add(p1z, p2z, xxt.get(2, 2));
	    	
	    	dfdx.add(p2x, p1x, xxt.get(0, 0));
	    	dfdx.add(p2x, p1y, xxt.get(0, 1));
	    	dfdx.add(p2x, p1z, xxt.get(0, 2));
	    	dfdx.add(p2y, p1x, xxt.get(1, 0));
	    	dfdx.add(p2y, p1y, xxt.get(1, 1));
	    	dfdx.add(p2y, p1z, xxt.get(1, 2));
	    	dfdx.add(p2z, p1x, xxt.get(2, 0));
	    	dfdx.add(p2z, p1y, xxt.get(2, 1));
	    	dfdx.add(p2z, p1z, xxt.get(2, 2));

    }   
 
    /**
     * Adds this springs damping contribution to the implicit damping matrix
     * @param dfdv
     */
    public void addDfdv( Matrix dfdv ) {
        // TODO: Objective 8, FINISH THIS CODE... necessary for backward euler integration
	    	DenseVector xhat = new DenseVector(3);
	    	DenseMatrix xxt = new DenseMatrix(3,3);
	    	
	    	// l = p1 - p2
	    	Vector3d l = new Vector3d();
		l.set(p1.p);
		l.sub(p2.p);
		l.normalize();
		
		// xhat = lhat 
		xhat.set(0,l.x);
		xhat.set(1,l.y);
		xhat.set(2,l.z);
		
		// xxt = x * xt
	    	xxt.rank1(xhat);
	    	// dfdv := - c * x * xt
	    	xxt.scale(-c);
	    	
	    	int p1x = 3 * p1.index;
	    	int p1y = 3 * p1.index + 1;
	    	int p1z = 3 * p1.index + 2;
	    	int p2x = 3 * p2.index;
	    	int p2y = 3 * p2.index + 1;
	    	int p2z = 3 * p2.index + 2;
	    	
	    	dfdv.add(p1x, p1x, xxt.get(0, 0));
	    	dfdv.add(p1x, p1y, xxt.get(0, 1));
	    	dfdv.add(p1x, p1z, xxt.get(0, 2));
	    	dfdv.add(p1y, p1x, xxt.get(1, 0));
	    	dfdv.add(p1y, p1y, xxt.get(1, 1));
	    	dfdv.add(p1y, p1z, xxt.get(1, 2));
	    	dfdv.add(p1z, p1x, xxt.get(2, 0));
	    	dfdv.add(p1z, p1y, xxt.get(2, 1));
	    	dfdv.add(p1z, p1z, xxt.get(2, 2));
	    	
	    	dfdv.add(p2x, p2x, xxt.get(0, 0));
	    	dfdv.add(p2x, p2y, xxt.get(0, 1));
	    	dfdv.add(p2x, p2z, xxt.get(0, 2));
	    	dfdv.add(p2y, p2x, xxt.get(1, 0));
	    	dfdv.add(p2y, p2y, xxt.get(1, 1));
	    	dfdv.add(p2y, p2z, xxt.get(1, 2));
	    	dfdv.add(p2z, p2x, xxt.get(2, 0));
	    	dfdv.add(p2z, p2y, xxt.get(2, 1));
	    	dfdv.add(p2z, p2z, xxt.get(2, 2));

	    	// dfdx_i = - dfdx_j
	    	xxt.scale(-1);
	    	dfdv.add(p1x, p2x, xxt.get(0, 0));
	    	dfdv.add(p1x, p2y, xxt.get(0, 1));
	    	dfdv.add(p1x, p2z, xxt.get(0, 2));
	    	dfdv.add(p1y, p2x, xxt.get(1, 0));
	    	dfdv.add(p1y, p2y, xxt.get(1, 1));
	    	dfdv.add(p1y, p2z, xxt.get(1, 2));
	    	dfdv.add(p1z, p2x, xxt.get(2, 0));
	    	dfdv.add(p1z, p2y, xxt.get(2, 1));
	    	dfdv.add(p1z, p2z, xxt.get(2, 2));
	    	
	    	dfdv.add(p2x, p1x, xxt.get(0, 0));
	    	dfdv.add(p2x, p1y, xxt.get(0, 1));
	    	dfdv.add(p2x, p1z, xxt.get(0, 2));
	    	dfdv.add(p2y, p1x, xxt.get(1, 0));
	    	dfdv.add(p2y, p1y, xxt.get(1, 1));
	    	dfdv.add(p2y, p1z, xxt.get(1, 2));
	    	dfdv.add(p2z, p1x, xxt.get(2, 0));
	    	dfdv.add(p2z, p1y, xxt.get(2, 1));
	    	dfdv.add(p2z, p1z, xxt.get(2, 2));
    } 
    
}
