package org.rcsb.assembly;

import javax.vecmath.Point3d;

/**
 * Base class for various vertex implementations, to work around
 * the lack of a KPartiteGraph implementation in JUNG2
 * 
 * Represents a vertex with a defined spatial location
 * @author spencer
 *
 */
public abstract class AtomVertex {
	protected Point3d position;
	public AtomVertex() {
		this.position = null;
	}
	public Point3d getPosition() {
		return position;
	}
	public void setPosition(Point3d pos) {
		this.position = pos;
	}
}
