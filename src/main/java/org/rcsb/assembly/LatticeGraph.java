package org.rcsb.assembly;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.contact.Pair;
import org.biojava.bio.structure.contact.StructureInterface;
import org.biojava.bio.structure.contact.StructureInterfaceList;
import org.biojava.bio.structure.gui.BiojavaJmol;
import org.biojava.bio.structure.xtal.CrystalBuilder;
import org.biojava.bio.structure.xtal.CrystalCell;
import org.biojava.bio.structure.xtal.CrystalTransform;
import org.biojava.bio.structure.xtal.SpaceGroup;
import org.biojava3.structure.quaternary.utils.Edge;
import org.biojava3.structure.quaternary.utils.SimpleGraph;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class LatticeGraph extends SimpleGraph<LatticeGraph.AtomVertex> {
	private static final Logger logger = LoggerFactory.getLogger(LatticeGraph.class);


	public static class AtomVertex {
		private String chainId;
		private Atom position;
		private int opId;
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


	public LatticeGraph(Structure struc) {

		// Get centroid for all chains in the asymmetric unit
		Map<String,Map<Integer,AtomVertex>> vertices = new HashMap<String, Map<Integer,AtomVertex>>();
		for(Chain c : struc.getChains()) {
			String chainId = c.getChainID();
			Atom[] ca = StructureTools.getAtomCAArray(c);
			Atom centroid = Calc.getCentroid(ca);

			AtomVertex vert = new AtomVertex(chainId,centroid,0);
			vertices.putIfAbsent(chainId, new HashMap<Integer,AtomVertex>());
			vertices.get(chainId).put(0, vert);

			this.addVertex(vert);
		}
		logger.info("Found "+vertices.size()+" chains in asymmetric unit");

		SpaceGroup space = struc.getCrystallographicInfo().getSpaceGroup();
		
		CrystalCell cell = struc.getCrystallographicInfo().getCrystalCell();
		// get all interfaces
		CrystalBuilder builder = new CrystalBuilder(struc);
		StructureInterfaceList interfaces = builder.getUniqueInterfaces();
		logger.info("Calculating ASA for "+interfaces.size()+" potential interfaces");
		interfaces.calcAsas();
		interfaces.removeInterfacesBelowArea();
		logger.info("Found "+interfaces.size()+" interfaces");

		// For each interface, add edge between centroids
		for(StructureInterface face : interfaces) {
			System.out.println(face);

			Pair<CrystalTransform> transforms = face.getTransforms();
			int idA = transforms.getFirst().getTransformId();
			int idB = transforms.getSecond().getTransformId();

			Pair<String> chainIds = face.getMoleculeIds();
			String chainA = chainIds.getFirst();
			String chainB = chainIds.getSecond();

			//			vertices.putIfAbsent(chainA, new HashMap<Integer,AtomVertex>());
			//			vertices.putIfAbsent(chainB, new HashMap<Integer,AtomVertex>());

			AtomVertex a = vertices.get(chainA).get(idA);
			AtomVertex b = vertices.get(chainB).get(idB);

			if(a == null) {
				Matrix4d crystalOp = transforms.getFirst().getMatTransform();
				Matrix4d realOp = cell.transfToOrthonormal(crystalOp);

				Atom pos = vertices.get(chainA).get(0).getPosition();
				Atom newPos = transformAtom(realOp, pos);

				a = new AtomVertex(chainA, newPos, idA);

				vertices.get(chainA).put(idA, a);
				addVertex(a);
			}
			if(b == null) {
				Matrix4d crystalOp = transforms.getSecond().getMatTransform();
				Matrix4d realOp = cell.transfToOrthonormal(crystalOp);

				Atom pos = vertices.get(chainB).get(0).getPosition();
				Atom newPos = transformAtom(realOp, pos);

				b = new AtomVertex(chainB, newPos, idB);

				vertices.get(chainB).put(idB, b);
				addVertex(b);
			}


			this.addEdge(a, b);
		}
	}


	private static Atom transformAtom(Matrix4d realOp, Atom pos) {
		Point3d posPt = new Point3d(pos.getCoords());
		realOp.transform(posPt); // The only meaningful line of this method
		double[] coords = new double[3];
		posPt.get(coords);
		Atom newPos = new AtomImpl();
		newPos.setCoords(coords);
		return newPos;
	}


	public String drawVertices() {
		StringBuilder str = new StringBuilder();
		for(AtomVertex vert : getVertices()) {
			Atom pos = vert.getPosition();
			//str.append(String.format("draw %s CIRCLE %f,%f,%f SCALE 1.0 DIAMETER 5.0; ", vert.getName(), pos.getX(),pos.getY(),pos.getZ() ));
			str.append(String.format("isosurface ID %s CENTER {%f,%f,%f} SPHERE 5.0;\n",
					vert.getName(), pos.getX(),pos.getY(),pos.getZ() ));
		}
		logger.info("JMOL:\n"+str.toString());
		return str.toString();
	}
	
	public String drawEdges() {
		StringBuilder str = new StringBuilder();
		for(Edge<AtomVertex> edge : getEdges()) {
			AtomVertex a = edge.getVertex1();
			AtomVertex b = edge.getVertex2();
			Atom posA = a.getPosition();
			Atom posB = b.getPosition();
			str.append(String.format("draw ID edge_%s_%s VECTOR {%f,%f,%f} {%f,%f,%f};\n",
					a.getName(),b.getName(),
					posA.getX(),posA.getY(),posA.getZ(),
					posB.getX()-posA.getX(),posB.getY()-posA.getY(),posB.getZ()-posA.getZ() ));
		}
		logger.info("JMOL:\n"+str.toString());
		return str.toString();
	}


	public static void main(String[] args) {
		String name;
		name = "1a99"; // See eppic-science #14
		name = "4MD1"; // rhodopsin, P63
		name = "1C8R"; // rhodopsin, P63, two trimer interfaces

		String filename = "/home/spencer/pdb/"+name.toLowerCase()+".pdb";

		try {
			Structure struc = StructureTools.getStructure(filename);

			LatticeGraph graph = new LatticeGraph(struc);

			BiojavaJmol jmol = new BiojavaJmol();
			//jmol.setStructure(struc);
			jmol.evalString(String.format("load \"%s\" {1 1 1};",filename));
			jmol.evalString("set unitcell {0 0 0};");
			//jmol.evalString("");
			// cartoon
			//jmol.evalString("hide null; select all;  spacefill off; wireframe off; backbone off; cartoon on;  select ligand; wireframe 0.16;spacefill 0.5; color cpk;  select *.FE; spacefill 0.7; color cpk ;  select *.CU; spacefill 0.7; color cpk ;  select *.ZN; spacefill 0.7; color cpk ;  select alls ON;");
			jmol.evalString("select all; spacefill off; wireframe off; backbone off; cartoon on; select none;");
			
			jmol.evalString(graph.drawVertices());
			jmol.evalString(graph.drawEdges());

		} catch (IOException e) {
			e.printStackTrace();
		} catch (StructureException e) {
			e.printStackTrace();
		}
	}
}
