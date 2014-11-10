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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseMultigraph;


public class LatticeGraph {
	private static final Logger logger = LoggerFactory.getLogger(LatticeGraph.class);


	protected static class ChainVertexKey {
		public int opId;
		public String chainId;
		public ChainVertexKey(String chainId,int opId) {
			this.opId = opId;
			this.chainId = chainId;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((chainId == null) ? 0 : chainId.hashCode());
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
			ChainVertexKey other = (ChainVertexKey) obj;
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
	
	private Graph<AtomVertex,InterfaceEdge> graph;
	private Map<ChainVertexKey, ChainVertex> chainNodes;


	public LatticeGraph(Structure struc) {
		graph = new UndirectedSparseMultigraph<AtomVertex, InterfaceEdge>();

		// Begin SciFi comments
		// SPACE OPS! Transform!
		Matrix4d[] spaceOps = struc.getCrystallographicInfo().getTransformationsOrthonormal();

		// GROUP OF SPACE OPS! Symmetry operations
		SpaceGroup space = struc.getCrystallographicInfo().getSpaceGroup();

		// CRYSTAL POWER CELL! Cell dimensions and angles
		CrystalCell cell = struc.getCrystallographicInfo().getCrystalCell();

		// Maps chainId and unit cell operator id to a vertex
		chainNodes = new HashMap<ChainVertexKey,ChainVertex>();


		// Generate vertices for unit cell
		for(Chain c : struc.getChains()) {
			// Calculate centroid position for this chain in the AU
			String chainId = c.getChainID();
			Atom[] ca = StructureTools.getAtomCAArray(c);
			Atom centroidAU = Calc.getCentroid(ca);

			for(int opId = 0; opId < spaceOps.length; opId++) {
				// Apply operator to centroid
				Atom centroid = (Atom)centroidAU.clone();
				Calc.transform(centroid, spaceOps[opId]);

				// Make sure it is inside the cell
				toUnitCell(centroid,cell);

				// Create new vertex & add to the graph
				ChainVertex vert = new ChainVertex(chainId,opId);
				vert.setPosition(centroid);

				chainNodes.put(new ChainVertexKey(chainId,opId), vert);
				graph.addVertex(vert);
			}

		}
		logger.info("Found "+graph.getVertexCount()+" chains in unit cell");
		
		
		// get all interfaces
		CrystalBuilder builder = new CrystalBuilder(struc);
		StructureInterfaceList interfaces = builder.getUniqueInterfaces();
		logger.info("Calculating ASA for "+interfaces.size()+" potential interfaces");
		interfaces.calcAsas(StructureInterfaceList.DEFAULT_ASA_SPHERE_POINTS/3, //fewer for performance
				Runtime.getRuntime().availableProcessors(),
				StructureInterfaceList.DEFAULT_MIN_COFACTOR_SIZE);
		interfaces.removeInterfacesBelowArea();
		logger.info("Found "+interfaces.size()+" interfaces");

		// For each interface, add edges for each asymm unit
		for(StructureInterface face : interfaces) {
			System.out.println(face);

			Pair<CrystalTransform> transforms = face.getTransforms();
			CrystalTransform transformA = transforms.getFirst();
			CrystalTransform transformB = transforms.getSecond();
//			int idA = transformA.getTransformId();
//			int idB = transformB.getTransformId();
			Matrix4d crystalTransformA = transformA.getMatTransform();
			Matrix4d crystalTransformB = transformB.getMatTransform();
			Matrix4d faceTransformA = cell.transfToOrthonormal(crystalTransformA);
			Matrix4d faceTransformB = cell.transfToOrthonormal(crystalTransformB);

			Pair<String> chainIds = face.getMoleculeIds();
			String chainA = chainIds.getFirst();
			String chainB = chainIds.getSecond();

			ChainVertex auA = chainNodes.get(new ChainVertexKey(chainA,0));
			ChainVertex auB = chainNodes.get(new ChainVertexKey(chainB,0));
			
			
			// Add edge for each asymmetric unit
			for(int opId = 0; opId < spaceOps.length; opId++) {
				// Calculate endpoints
				// First transform each centroid according to the spaceOp (cached in the vertices)
				Atom startPosA = chainNodes.get(new ChainVertexKey(chainA,opId)).getPosition();
				Atom startPosB = chainNodes.get(new ChainVertexKey(chainB,opId)).getPosition();
				// Then transform according to the interface
				Atom endPosA = transformAtom(faceTransformA, startPosA);
				Atom endPosB = transformAtom(faceTransformB, startPosB);
				// Return to the Unit cell
				Atom ucPosA = (Atom) endPosA.clone();
				Atom ucPosB = (Atom) endPosB.clone();
				toUnitCell(ucPosA, cell);
				toUnitCell(ucPosB, cell);
				
				// Determine which AU the partners belong to
				ChainVertex vertA = findVertex(ucPosA);
				ChainVertex vertB = findVertex(ucPosB);
				
				//TODO Add edge
			}
//			AtomVertex a = vertices.get(chainA).get(idA);
//			AtomVertex b = vertices.get(chainB).get(idB);
//
//			if(a == null) {
//				Matrix4d crystalOp = transforms.getFirst().getMatTransform();
//				Matrix4d realOp = cell.transfToOrthonormal(crystalOp);
//
//				Atom pos = vertices.get(chainA).get(0).getPosition();
//				Atom newPos = transformAtom(realOp, pos);
//
//				a = new AtomVertex(chainA, idA);
//				a.setPosition(newPos);
//
//				vertices.get(chainA).put(idA, a);
//				graph.addVertex(a);
//			}
//			if(b == null) {
//				Matrix4d crystalOp = transforms.getSecond().getMatTransform();
//				Matrix4d realOp = cell.transfToOrthonormal(crystalOp);
//
//				Atom pos = vertices.get(chainB).get(0).getPosition();
//				Atom newPos = transformAtom(realOp, pos);
//
//				b = new AtomVertex(chainB, idB);
//				b.setPosition(newPos);
//
//				vertices.get(chainB).put(idB, b);
//				graph.addVertex(b);
//			}
//
//			InterfaceEdge edge = new InterfaceEdge(face.getId());
//			edge.addSegment(start, end);
//			
//			graph.addEdge(edge,a, b);
		}
	}

	/**
	 * Finds a Vertex that corresponds to the specified atom. 
	 * 
	 * Returns null if no vertex is found within a small margin of error
	 * @param vertices
	 * @param atom
	 * @return
	 */
	private ChainVertex findVertex(Atom atom) {
		final double tol = 1e-12;
		
		for(ChainVertex vert : chainNodes.values()) {
			if(vert.getPosition() == null) {
				continue;
			}
			double distSq = Calc.getDistanceFast(vert.getPosition(), atom);
			if(distSq < tol) {
				return vert;
			}
		}
		return null;
	}

	/**
	 * Transforms an Atom according to a matrix.
	 * 
	 * @param realOp A 4D rotation matrix
	 * @param pos The Atom to rotate (unmodified)
	 * @return A new Atom with the transformed coordinates
	 * @see Matrix4d#transform(Point3d)
	 */
	private static Atom transformAtom(Matrix4d realOp, Atom pos) {
		Point3d posPt = new Point3d(pos.getCoords());
		realOp.transform(posPt); // The only meaningful line of this method
		double[] coords = new double[3];
		posPt.get(coords);
		Atom newPos = new AtomImpl();
		newPos.setCoords(coords);
		return newPos;
	}
	
	/**
	 * Modify an Atom so that it resides within the origin unit cell
	 * @param atom The atom to be modified
	 * @param cell Definition of the unit cell
	 * @see CrystalCell#transfToOriginCell(javax.vecmath.Tuple3d)
	 */
	private static void toUnitCell(Atom atom, CrystalCell cell) {
		double[] coords = atom.getCoords();
		Point3d pt = new Point3d(coords);
		cell.transfToOriginCell(pt);
		atom.setX(pt.getX());
		atom.setY(pt.getY());
		atom.setZ(pt.getZ());
	}


	public String drawVertices() {
		StringBuilder str = new StringBuilder();
		for(AtomVertex vertex : graph.getVertices()) {
			if( vertex instanceof ChainVertex) {
				ChainVertex vert = (ChainVertex) vertex;
				Atom pos = vert.getPosition();
				//str.append(String.format("draw %s CIRCLE %f,%f,%f SCALE 1.0 DIAMETER 5.0; ", vert.getName(), pos.getX(),pos.getY(),pos.getZ() ));
				str.append(String.format("isosurface ID %s CENTER {%f,%f,%f} SPHERE 5.0;\n",
						vert, pos.getX(),pos.getY(),pos.getZ() ));
			} else {
				throw new IllegalStateException("Unrecognized vertex class "+vertex.getClass().toString());
			}
		}
		logger.info("JMOL:\n"+str.toString());
		return str.toString();
	}

	public String drawEdges() {
		StringBuilder str = new StringBuilder();
		for(InterfaceEdge edge : graph.getEdges()) {
			edu.uci.ics.jung.graph.util.Pair<AtomVertex> edgePair = graph.getEndpoints(edge);
			AtomVertex a = edgePair.getFirst();
			AtomVertex b = edgePair.getSecond();
			Atom posA = a.getPosition();
			Atom posB = b.getPosition();
			str.append(String.format("draw ID edge_%s_%s VECTOR {%f,%f,%f} {%f,%f,%f};\n",
					a,b,
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
