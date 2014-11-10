package org.rcsb.assembly;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

import edu.uci.ics.jung.graph.util.Pair;


class InterfaceEdge {
	// annotation data
	private int interfaceId;
	private List<Pair<Point3d>> segments;
	
	public InterfaceEdge(int interfaceId) {
		this.interfaceId = interfaceId;
		segments = new ArrayList<Pair<Point3d>>(2);
		}
	
	public void addSegment(Pair<Point3d> segment) {
		segments.add(segment);
	}
	public void addSegment(Point3d start, Point3d end) {
		addSegment(new Pair<Point3d>(start,end));
	}
	public List<Pair<Point3d>> getSegments() {
		return segments;
	}
	public int getInterfaceId() {
		return interfaceId;
	}
	
	@Override
	public String toString() {
		return String.format("-%d-",interfaceId);
	}
}