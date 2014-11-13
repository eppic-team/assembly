package org.rcsb.assembly;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JFrame;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Point3i;

import org.biojava.bio.structure.Atom;
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
	// Maps chainId and unit cell operator id to a vertex
	private Map<ChainVertexKey, ChainVertex> chainNodes;


	public LatticeGraph(Structure struc) {
		graph = new UndirectedSparseMultigraph<AtomVertex, InterfaceEdge>();
		chainNodes = new HashMap<ChainVertexKey,ChainVertex>();

		
		// Begin SciFi comments
		// SPACE OPS! Transform!
		Matrix4d[] spaceOps = struc.getCrystallographicInfo().getTransformationsOrthonormal();

		// GROUP OF SPACE OPS! Symmetry operations
		SpaceGroup spaceGroup = struc.getCrystallographicInfo().getSpaceGroup();

		// CRYSTAL POWER CELL! Cell dimensions and angles
		CrystalCell cell = struc.getCrystallographicInfo().getCrystalCell();
		
		

		// Generate vertices for unit cell
		initChainVertices(struc, spaceOps, cell);
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
			if(face.getId() != 2) continue;
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
			
			Point3d startPosA = auA.getPosition();
			Point3d startPosB = auB.getPosition();

			// transform according to the interface
			Point3d endPosA = new Point3d(startPosA);
			Point3d endPosB = new Point3d(startPosB);
			faceTransformA.transform(endPosA);
			faceTransformB.transform(endPosB);

			
			// Add edge for each asymmetric unit
			for(int opId = 0; opId < spaceOps.length; opId++) {
				//if(opId != 5) continue;
				
				// transform according to the spaceop
				Point3d posA = new Point3d(endPosA);
				Point3d posB = new Point3d(endPosB);
				spaceOps[opId].transform(posA);
				spaceOps[opId].transform(posB);

//				
//				// Calculate endpoints
//				// First transform each centroid according to the spaceOp (cached in the vertices)
//				Atom startPosA = chainNodes.get(new ChainVertexKey(chainA,opId)).getPosition();
//				Atom startPosB = chainNodes.get(new ChainVertexKey(chainB,opId)).getPosition();
//				// Then transform according to the interface
//				Atom endPosA = transformAtom(faceTransformA, startPosA);
//				Atom endPosB = transformAtom(faceTransformB, startPosB);
				// Return to the Unit cell
				Point3d ucPosA = new Point3d(posA);
				Point3d ucPosB = new Point3d(posB);
				cell.transfToOriginCell(ucPosA);
				cell.transfToOriginCell(ucPosB);
				logger.info(String.format("Interface %d AU %d: %s;\t%s\t->\t%s;\t%s%n",face.getId(),opId,
						round(posA,2),round(posB,2),round(ucPosA,2),round(ucPosB,2) ));

				
				// Determine which AU the partners belong to
				ChainVertex vertA = findVertex(ucPosA);
				ChainVertex vertB = findVertex(ucPosB);
				
				//TODO remove duplicates
				// Only draw edges where cell(A) <= cell(B)
				Point3i cellA = cell.getCellIndices(endPosA);
				Point3i cellB = cell.getCellIndices(endPosB);
				if( cellA.x > cellB.x || 
						cellA.y > cellB.y || 
						cellA.z > cellB.z )
				{
					//continue;
				}
				// If equal, draw edges where au(A) < au(B)
				else if( cellA.equals(cellB) && transformA.getTransformId() > transformB.getTransformId()) {
					//continue;
				}
				// Use chainId as final tie breaker
				else if(transformA.getTransformId() ==transformB.getTransformId() &&
						chainA.compareTo(chainB) > 0) {
					//continue;
				}

				
				// Create interface vertex at the midpoint
				InterfaceVertex ivert = new InterfaceVertex(opId,face.getId());
				Point3d mid = new Point3d();
				mid.add(posA, posB);
				mid.scale(0.5);
				Point3d ucMid = new Point3d(mid);
				cell.transfToOriginCell(ucMid);
				ivert.setPosition(ucMid);
				graph.addVertex(ivert);
				
				InterfaceEdge edgeA = new InterfaceEdge(face.getId());
				InterfaceEdge edgeB = new InterfaceEdge(face.getId());
				edgeA.setColor("green");
				edgeB.setColor("blue");
				
				//TODO Set segments for wrapped edges
				
				
				graph.addEdge(edgeA, vertA, ivert);
				graph.addEdge(edgeB, vertB, ivert);
			}
			//break;
			
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

	private Point3d round(Point3d p, int places) {
		double placeMult = Math.pow(10, places);
		double x = Math.round(p.x*placeMult)/placeMult;
		double y = Math.round(p.y*placeMult)/placeMult;
		double z = Math.round(p.z*placeMult)/placeMult;
		return new Point3d( x,y,z );
	}

	/**
	 * Initialize the ChainVertex nodes of the graph
	 * @param struc
	 * @param spaceOps
	 * @param cell
	 */
	private void initChainVertices(Structure struc, Matrix4d[] spaceOps,
			CrystalCell cell) {
		for(Chain c : struc.getChains()) {
			// Calculate centroid position for this chain in the AU
			String chainId = c.getChainID();
			Atom[] ca = StructureTools.getAtomCAArray(c);
			Point3d centroidAU = new Point3d(Calc.getCentroid(ca).getCoords());

			for(int opId = 0; opId < spaceOps.length; opId++) {
				// Apply operator to centroid
				Point3d centroid = new Point3d(centroidAU);
				spaceOps[opId].transform(centroid);

				// Make sure it is inside the cell
				cell.transfToOriginCell(centroid);

				// Create new vertex & add to the graph
				ChainVertex vert = new ChainVertex(chainId,opId);
				vert.setPosition(centroid);

				chainNodes.put(new ChainVertexKey(chainId,opId), vert);
				graph.addVertex(vert);
			}

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
	private ChainVertex findVertex(Point3d atom) {
		final double tol = 1e-12;
		
		for(ChainVertex vert : chainNodes.values()) {
			if(vert.getPosition() == null) {
				continue;
			}
			double distSq = vert.getPosition().distanceSquared(atom);
			if(distSq < tol) {
				return vert;
			}
		}
		return null;
	}

	public String drawVertices() {
		StringBuilder str = new StringBuilder();
		for(AtomVertex vertex : graph.getVertices()) {
			if( vertex instanceof ChainVertex) {
				ChainVertex vert = (ChainVertex) vertex;
				Point3d pos = vert.getPosition();
				str.append(String.format("isosurface ID chain%s CENTER {%f,%f,%f} SPHERE 5.0 COLOR blue;\n",
						vert, pos.x,pos.y,pos.z ));
				str.append(String.format("set echo ID echoChain%s {%f,%f,%f}; color echo blue; echo \"  %s\";\n",
						vert, pos.x,pos.y,pos.z,vert));
			} else if( vertex instanceof InterfaceVertex ) {
				InterfaceVertex vert = (InterfaceVertex) vertex;
				
				Point3d pos = vert.getPosition();
				
				str.append(String.format("draw ID interface%s CIRCLE {%f,%f,%f} DIAMETER 5.0;\n",vert,
						pos.x,pos.y,pos.z ));
				str.append(String.format("set echo ID echoInterface%s {%f,%f,%f}; color echo green;  echo %s;\n",
						vert, pos.x,pos.y,pos.z,vert.getInterfaceId()));
			} else {
				throw new IllegalStateException("Unrecognized vertex class "+vertex.getClass().toString());
			}
		}
		logger.trace("JMOL:\n"+str.toString());
		return str.toString();
	}

	public String drawEdges() {
		StringBuilder str = new StringBuilder();
		for(InterfaceEdge edge : graph.getEdges()) {
			edu.uci.ics.jung.graph.util.Pair<AtomVertex> edgePair = graph.getEndpoints(edge);
			AtomVertex a = edgePair.getFirst();
			AtomVertex b = edgePair.getSecond();
			Point3d posA = a.getPosition();
			Point3d posB = b.getPosition();
			
			double xjitter = (Math.random()*2-1.0) * 1.0;
			double yjitter = (Math.random()*2-1.0) * 1.0;
			double zjitter = (Math.random()*2-1.0) * 1.0;
			str.append(String.format("draw ID edge_%s_%s VECTOR {%f,%f,%f} {%f,%f,%f} COLOR %s;\n",
					a,b,
					posA.x+xjitter,posA.y+yjitter,posA.z+zjitter,
					posB.x-posA.x,posB.y-posA.y,posB.z-posA.z,
					edge.getColor() == null ? "yellow" : edge.getColor() ));
		}
		logger.trace("JMOL:\n"+str.toString());
		return str.toString();
	}


	public static void main(String[] args) {
		String name;
		name = "1a99"; // See eppic-science #14
		name = "4MD1"; // rhodopsin, P63
		name = "1C8R"; // rhodopsin, P63, two trimer interfaces

		String filename = System.getProperty("user.home")+"/pdb/"+name.toLowerCase()+".pdb";

		try {
			Structure struc = StructureTools.getStructure(filename);

			LatticeGraph graph = new LatticeGraph(struc);

			BiojavaJmol jmol = new BiojavaJmol();
			//jmol.setStructure(struc);
			jmol.evalString(String.format("load \"%s\" {1 1 1};",filename));
			jmol.evalString("set unitcell {0 0 0};");
			//jmol.evalString("");
			jmol.evalString(graph.drawVertices());
			jmol.evalString(graph.drawEdges());
			// cartoon
			//jmol.evalString("hide null; select all;  spacefill off; wireframe off; backbone off; cartoon on;  select ligand; wireframe 0.16;spacefill 0.5; color cpk;  select *.FE; spacefill 0.7; color cpk ;  select *.CU; spacefill 0.7; color cpk ;  select *.ZN; spacefill 0.7; color cpk ;  select alls ON;");
			jmol.evalString("select all; spacefill off; wireframe off; backbone off; cartoon off; select none;");
			jmol.evalString("set axes molecular;");

			jmol.getFrame().setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		} catch (IOException e) {
			e.printStackTrace();
		} catch (StructureException e) {
			e.printStackTrace();
		}
	}
}
