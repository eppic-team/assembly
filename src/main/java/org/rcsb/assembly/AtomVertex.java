package org.rcsb.assembly;

import org.biojava.bio.structure.Atom;

/**
 * Base class for various vertex implementations, to work around
 * the lack of a KPartiteGraph implementation in JUNG2
 * 
 * Represents a vertex with a defined spatial location
 * @author spencer
 *
 */
public abstract class AtomVertex {
	protected Atom position;
	public AtomVertex() {
		this.position = null;
	}
	public Atom getPosition() {
		return position;
	}
	public void setPosition(Atom pos) {
		this.position = pos;
	}
}
