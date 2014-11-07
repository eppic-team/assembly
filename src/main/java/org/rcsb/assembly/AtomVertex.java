package org.rcsb.assembly;

import org.biojava.bio.structure.Atom;

public class AtomVertex {
	// Primary Key:
	private int opId; // operator to generate this position within the unit cell
	private String chainId;
	
	// Metadata
	private Atom position;
	private int entity; //TODO
	
	public AtomVertex(String chainId, Atom position, int opId) {
		super();
		this.chainId = chainId;
		this.position = position;
		this.opId = opId;
	}
	public String getName() {
		return chainId+opId;
	}
	public String getChainId() {
		return chainId;
	}
	public Atom getPosition() {
		return position;
	}
	public int getOpId() {
		return opId;
	}
	
	@Override
	public String toString() {
		return getName();
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chainId == null) ? 0 : chainId.hashCode());
		result = prime * result + opId;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		AtomVertex other = (AtomVertex) obj;
		if (chainId == null) {
			if (other.chainId != null)
				return false;
		} else if (!chainId.equals(other.chainId))
			return false;
		if (opId != other.opId)
			return false;
		return true;
	}

}