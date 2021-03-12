// NAME: Zhiguo(Frank) Zhang
// ID: 260550226
package comp559.project;

import java.util.LinkedList;
import java.util.List;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * Implementation of a robust collision detection.
 * @author kry
 */
public class RobustCCD {

    double restitution;
    
    double H;
    
    int MAX_ITERATION;
    
    /**
     * Creates the new continuous collision detection and response object
     */
    public RobustCCD() {
        // do nothing
    }
    
    /**
     * Checks all collisions in interval t to t+h
     * @param h
     * @param system 
     * @return true if all collisions resolved
     */
    public boolean check( double h, ParticleSystem system ) {        

        // For each particle-edge pair, find the roots for when the three particles are
        // co-linear, and then pick the first root on (0,h] which corresponds to an 
    		// actual collision.  Compute a collision response.  That is, compute an appropriate
    		// collision normal, compute the impulse, and then apply the impulse to the associated
    		// particles.  Be sure to deal with pinning constraints!  Repeat until all collisions
    		// are resolved and it is safe to advance time
       
     	// You probably want to write other methods to help with the job, and
     	// call them here.  Or alternatively you can write one large and ugly
     	// monolithic function here.
  
        // TODO: OBJECTIVE 1 continuous collision detection    	
    		// TODO: OBJECTIVE 2 compute collision impulses
    		// TODO: OBJECTIVE 3 iterative to resolve all collisions
		return iterativeCheck(h, system, 0);
    }
    private List<Particle> allParticles = null;
    private List<Spring> allSprings = null;
    private List<Triangle> allTriangles = null;
    /**
     * Iteratively resolve all collisions, give up if MAX_ITERATION is reached
     * @param h Time step
     * @param system Particle-spring system
     * @param iter current iteration
     * @return true if all collisions are resolved, false if give up
     */
    private boolean iterativeCheck(double h, ParticleSystem system, int iter) {
		if (iter > MAX_ITERATION) return false; // give up
		
	    	boolean collision = false;
	    	
	    	allParticles = new LinkedList<Particle>(system.particles);
	    	allParticles.addAll(system.objParticles);
	    	
	    	allSprings = new LinkedList<Spring>(system.springs);
	    	allSprings.addAll(system.objSprings);
	    	
	    	allTriangles = new LinkedList<Triangle>(system.triangles);
	    	allTriangles.addAll(system.objTriangles);
	    	
	    	for(Particle p: allParticles) {
	    		for(Triangle t: allTriangles) {
	    			if(p != t.p1 && p != t.p2 && p != t.p3) {
					Point3d p1 = t.p1.p;
					Point3d p2 = t.p2.p;
					Point3d p3 = t.p3.p;
					Vector3d v1 = new Vector3d();
					Vector3d v2 = new Vector3d();
					Vector3d normal = new Vector3d();

					v1.sub(p2,p1);
					v2.sub(p3,p1);
					normal.cross(v1,v2);
					normal.normalize();

					Vector3d v = new Vector3d();
					v.sub(p.p,p1);

					double dot = normal.dot(v);
					if(Math.abs(dot) < 1e-6) {
						double l = v.dot(normal);
						normal.scale(l);

						Point3d project = new Point3d();
						project.sub(p.p, normal);

						Vector3d v0 = new Vector3d();
						v0.sub(project,p1);

						Vector3d normalA = new Vector3d();
						Vector3d normalB = new Vector3d();
						Vector3d normalC = new Vector3d();

						Vector3d v3 = new Vector3d();
						Vector3d v4 = new Vector3d();
						Vector3d v5 = new Vector3d();
						Vector3d v6 = new Vector3d();

						v3.sub(p3,p2);
						v4.sub(project,p2);
						v5.sub(p1,p3);
						v6.sub(project,p3);

						normalA.cross(v3, v4);
						normalB.cross(v5, v6);
						normalC.cross(v1, v0);

						normal.cross(v1,v2);

						double a = normal.dot(normalA)/normal.lengthSquared();
						double b = normal.dot(normalB)/normal.lengthSquared();
						double c = normal.dot(normalC)/normal.lengthSquared();

						normal.normalize();

						if(0 < a && a < 1 && 0 < b && b < 1 && 0 < c && c < 1) {
							Vector3d temp = new Vector3d();
							temp.sub(project,p.p);

							if(normal.dot(p.v) < 0) normal.negate();


							Vector3d vRel = new Vector3d();
							Vector3d vTriangle = new Vector3d(a*t.p1.v.x + b*t.p2.v.x + c*t.p3.v.x, a*t.p1.v.y + b*t.p2.v.y + c*t.p3.v.y, a*t.p1.v.z + b*t.p2.v.z + c*t.p3.v.z);

							vRel.sub(p.v, vTriangle);
							double vrel_neg = vRel.dot(normal);
							
							double i =  ( 1 + restitution ) * vrel_neg / (a * a / t.p1.mass + b * b / t.p2.mass+ c * c / t.p3.mass+ 1.0/p.mass);
									
							if (!p.pinned) p.v.add(new Vector3d(-i*normal.x/p.mass, -i*normal.y/p.mass, -i*normal.z/p.mass));
							if (!t.p1.pinned) t.p1.v.add(new Vector3d(i*a*normal.x/t.p1.mass, i*a*normal.y/t.p1.mass, i*a*normal.z/t.p1.mass));
							if (!t.p2.pinned) t.p2.v.add(new Vector3d(i*b*normal.x/t.p2.mass, i*b*normal.y/t.p2.mass, i*b*normal.z/t.p2.mass));
							if (!t.p3.pinned) t.p3.v.add(new Vector3d(i*c*normal.x/t.p3.mass, i*c*normal.y/t.p3.mass, i*c*normal.z/t.p3.mass));
							collision = true;
						}
					}
	    			}
			}
	    	}
	    	
	    	// edge vs edge test
		for (Spring s1 : allSprings) {
			for (Spring s2 : allSprings) {
				if ( s1 != s2 ) {
					// check if two springs intersect
					// http://paulbourke.net/geometry/pointlineplane/
					double EPS = 1e-6;
					Point3d p1, p2, p3, p4, pa, pb, p13, p43, p21;
					double mua, mub, d1343, d4321, d1321, d4343, d2121, numer, denom;
					p1 = s1.p1.p;
					p2 = s1.p2.p;
					p3 = s2.p1.p;
					p4 = s2.p2.p;
					p13 = new Point3d();
					p13.sub(p1, p3);
					p43 = new Point3d();
					p43.sub(p4, p3);
					if (Math.abs(p43.x) < EPS && Math.abs(p43.y) < EPS && Math.abs(p43.z) < EPS) continue;
					p21 = new Point3d();
					p21.sub(p2, p1);
					if (Math.abs(p21.x) < EPS && Math.abs(p21.y) < EPS && Math.abs(p21.z) < EPS) continue;
					d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
					d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
					d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
					d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
					d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;

					denom = d2121 * d4343 - d4321 * d4321;
					if (Math.abs(denom) < EPS) continue;
					numer = d1343 * d4321 - d1321 * d4343;

					mua = numer / denom;
					mub = (d1343 + d4321 * (mua)) / d4343;
					if (mua < 0 || mua > 1 || mub < 0 || mub > 1) continue;
					pa = new Point3d();
					pa.scaleAdd(mua, p21, p1);
					pb = new Point3d();
					pb.scaleAdd(mub, p43, p3);

					// TODO: finish collision
					Vector3d pab = new Vector3d();
					pab.sub(pb, pa);
					if (pab.length() < EPS) {
						Vector3d vRel = new Vector3d();
						Vector3d v1 = new Vector3d(mua * s1.p2.v.x + (1-mua) * s1.p1.v.x, mua * s1.p2.v.y + (1-mua) * s1.p1.v.y, mua * s1.p2.v.z + (1-mua) * s1.p1.v.z);
						Vector3d v2 = new Vector3d(mub * s2.p2.v.x + (1-mub) * s2.p1.v.x, mub * s2.p2.v.y + (1-mub) * s2.p1.v.y, mub * s2.p2.v.z + (1-mub) * s2.p1.v.z);

						vRel.sub(v1, v2);
						Vector3d normal = new Vector3d(pab);
						normal.normalize();
						double vrel_neg = vRel.dot(normal);
						double i =  ( 1 + restitution ) * vrel_neg / (mua*mua/s1.p2.mass + (1-mua)*(1-mua)/s1.p1.mass + mub*mub/s2.p2.mass + (1-mub)*(1-mub)/s2.p1.mass);

						if (!s1.p1.pinned) s1.p1.v.add(new Vector3d(-(1-mua)*i*normal.x/s1.p1.mass,-(1-mua)*i*normal.y/s1.p1.mass,-(1-mua)*i*normal.z/s1.p1.mass));
						if (!s1.p2.pinned) s1.p2.v.add(new Vector3d(  -(mua)*i*normal.x/s1.p2.mass,  -(mua)*i*normal.y/s1.p2.mass,  -(mua)*i*normal.z/s1.p2.mass));
						if (!s2.p1.pinned) s2.p1.v.add(new Vector3d( (1-mub)*i*normal.x/s2.p1.mass, (1-mub)*i*normal.y/s2.p1.mass, (1-mub)*i*normal.z/s2.p1.mass));
						if (!s2.p2.pinned) s2.p2.v.add(new Vector3d(   (mub)*i*normal.x/s2.p2.mass,   (mub)*i*normal.y/s2.p2.mass,   (mub)*i*normal.z/s2.p2.mass));
						collision = true;
					}
				}
			}
		}
		
	    	if (collision) return iterativeCheck(h, system, iter+1); // Might be more collision!
	    	else return true; // No more collision found!
    }
}
