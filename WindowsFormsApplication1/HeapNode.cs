using System;
using System.Collections;
namespace TSP
{
	// data structure used in both Array and Heap implementations of PriorityQueue. Holds a current distance, the index 
	// of the node, and the index of the previous point. for this node.
	public class HeapNode
	{
		public ArrayList CurrentPath;
		// Space Compleixty of O(n^2)
		public double[,] CostArray;
		public double LowerBound;
		public int LastCity;
		public int NumberOfNodes;

		// Space Complexity of O(n^2)
		// Each HeapNode requires O(n^2) because each node has a cost matrix that is n by n in size
		public HeapNode(double[,] costArray, ArrayList currentPath, double lowerBound)
		{
			this.CostArray = costArray;
			this.CurrentPath = currentPath;
			this.LastCity = (int) currentPath[currentPath.Count - 1];
			this.LowerBound = lowerBound;
			this.NumberOfNodes = currentPath.Count;
		}

		// Helper function to get the cost from the last node visited in this partial solution
		// to a given next city index
		public double costToCity(int cityIndex) {
			return this.CostArray[this.LastCity, cityIndex];
		}
	}
}
