package org.rcsb.assembly;

import java.awt.Dimension;

import javax.swing.JFrame;

import edu.uci.ics.jung.algorithms.layout.CircleLayout;
import edu.uci.ics.jung.algorithms.layout.Layout;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;
import edu.uci.ics.jung.visualization.BasicVisualizationServer;
import edu.uci.ics.jung.visualization.decorators.ToStringLabeller;

public class GraphTest {

	public static void main(String[] args) {
		UndirectedSparseGraph<ChainVertex, InterfaceEdge> graph = new UndirectedSparseGraph<ChainVertex, InterfaceEdge>();

		ChainVertex v1 = new ChainVertex("A", 1);
		ChainVertex v2 = new ChainVertex("B", 1);
		ChainVertex v3 = new ChainVertex("C", 1);
		graph.addVertex(v1);
		graph.addVertex(v2);
		graph.addVertex(v3);

		ChainVertex fake1 = new ChainVertex("A", 1);
		ChainVertex fake2 = new ChainVertex("B", 1);

		InterfaceEdge e1 = new InterfaceEdge(1);
		InterfaceEdge e2 = new InterfaceEdge(2);

		graph.addEdge(e1, fake1,fake2);
		graph.addEdge(e2, v2, v3);

		Pair<ChainVertex> ends = graph.getEndpoints(e1);
		System.out.format("v1 %s first %s fake1%n", v1 == ends.getFirst() ? "==" : "!=", ends.getFirst() == fake1 ? "==" : "!=");
		System.out.format("v2 %s second %s fake2%n", v2 == ends.getSecond() ? "==" : "!=", ends.getSecond() == fake2 ? "==" : "!=");
		ends = graph.getEndpoints(e2);
		System.out.format("v2 %s first %s fake2%n", v2 == ends.getFirst() ? "==" : "!=", ends.getFirst() == fake2 ? "==" : "!=");

		System.out.format("Graph has %d vertices and %d edges.%n",graph.getVertexCount(),graph.getEdgeCount());

		// Visualize


		Layout<ChainVertex, InterfaceEdge> layout = new CircleLayout<ChainVertex, InterfaceEdge>(graph);
		layout.setSize(new Dimension(300,300)); // sets the initial size of the space
		
		// The BasicVisualizationServer<V,E> is parameterized by the edge types
		BasicVisualizationServer<ChainVertex, InterfaceEdge> vv =
				new BasicVisualizationServer<ChainVertex, InterfaceEdge>(layout);

		vv.getRenderContext().setVertexLabelTransformer(new ToStringLabeller<ChainVertex>());
		vv.getRenderContext().setEdgeLabelTransformer(new ToStringLabeller<InterfaceEdge>());

		vv.setPreferredSize(new Dimension(350,350)); //Sets the viewing area size
		JFrame frame = new JFrame("Simple Graph View");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.getContentPane().add(vv);
		frame.pack();
		frame.setVisible(true);
	}
}
