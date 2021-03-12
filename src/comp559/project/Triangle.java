// NAME: Zhiguo(Frank) Zhang
// ID: 260550226
package comp559.project;

import javax.vecmath.Vector3d;
public class Triangle {
	
	Particle p1 = null;
	Particle p2 = null;
	Particle p3 = null;
	
	Vector3d n = new Vector3d();
	public Triangle(Particle p1, Particle p2, Particle p3) {
		this.p1 = p1;
        this.p2 = p2;
        this.p3 = p3;
	}
	
	public Vector3d normal() {
		Vector3d ab = new Vector3d();
		Vector3d ac = new Vector3d();
		ab.sub(p2.p, p1.p);
		ac.sub(p3.p, p1.p);
		this.n.cross(ab, ac);
		this.n.normalize();
		return this.n;
	}
	public double distance(double x, double y, double z) {
		normal();
		double d = - n.dot(new Vector3d(this.p1.p));
		Vector3d p = new Vector3d(x,y,z);
		return (p.dot(n) + d) / n.length();
	}
}
