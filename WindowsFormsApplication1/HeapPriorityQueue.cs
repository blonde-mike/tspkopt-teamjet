using System;
using System.Drawing;
using System.Collections.Generic;

namespace TSP
{
	public class HeapPriorityQueue
	{
		// nodes is the heapified array. The positions of each node swaps on each insert, decreaseKey, or deleteMin
		List<HeapNode> nodes;
		public int count = 0;
		public int maxCountEver = 0;
		public HeapPriorityQueue()
		{
			nodes = new List<HeapNode>();
		}

		// Time Complexity: O(log n) decreaseKey is the key of keeping a binary heap balanced. It updates the position 
		// of a node at a given key by bubbling it up through the tree. It compares the node's distance value with the
		// value of its parent, and if the node is less than it, then we swap them. We run this algorithm until node has
		// become the top of the min binary heap (position 0), or until its parent's distance value is less than its own.
		// The worst case of this function is if the updated node at a given key is at the bottom of the heap, and bubbles
		// up to the very top of the heap. Because our heap is a balanced binary tree, the max height is log n (n 
		// being the number of nodes in the tree) thus the maximum this function will iterate is O(log n).
		private void decreaseKey(int key, HeapNode pointNode)
		{
			int parent;
			int newIndex = key;
			while (newIndex > 0)
			{
				parent = (int)(Math.Ceiling(Convert.ToDouble(newIndex) / 2)) - 1;

				if ((pointNode.NumberOfNodes > nodes[parent].NumberOfNodes) ||  (pointNode.NumberOfNodes == nodes[parent].NumberOfNodes && pointNode.LowerBound < nodes[parent].LowerBound))
				{
					//Console.WriteLine("swapping");
					swap(newIndex, parent);
					newIndex = parent;
				}
				else
				{
					break;
				}

			}
			// Console.WriteLine("done decreasingKey");
		}

		// Time Complexity: O(log n). deleteMin removes the top node in O(1), but in order to keep the binary heap 
		// balanced, it moves the last element of the heap to the first position then bubbles it down through the heap.
		// Because our heap is a balanced binary tree, the max height is log n (n being the number of nodes in the tree),
		// thus the maximum this function will iterate is O(log n).
		public HeapNode deleteMin()
		{
			HeapNode returner = nodes[0];
			HeapNode newLast = nodes[count - 1];
			nodes[0] = newLast;
			count -= 1;

			int current = 0;
			while (current < count)
			{
				int leftChild = current * 2 + 1;

				// if leftChild >= count, then neither child exists
				if (leftChild >= count)
				{
					break;
				}

				int leftCount = nodes[leftChild].NumberOfNodes;
				double leftDist = nodes[leftChild].LowerBound;
				int rightCount = int.MaxValue;
				double rightDist = double.MaxValue;

				int rightChild = leftChild + 1;

				// if rightChild is valid
				if (rightChild < count)
				{
					rightDist = nodes[rightChild].LowerBound;
					rightCount = nodes[rightChild].NumberOfNodes;
				}

				if ( (leftCount > rightCount) ||  (leftCount == rightCount && leftDist < rightDist))
				{
					if ( (leftCount > newLast.NumberOfNodes) || (leftCount == newLast.NumberOfNodes && leftDist < newLast.LowerBound) )
					{
						// left child is best
						swap(leftChild, current);
						current = leftChild;
						continue;
					}
				}
				else
				{
					if ( (rightCount > newLast.NumberOfNodes) || (rightCount == newLast.NumberOfNodes && rightDist < newLast.LowerBound) )
					{
						swap(rightChild, current);
						current = rightChild;
						continue;
					}
				}

				// if this line is reached, then our newLast was either greater than its smaller child.
				break;
			}

			return returner;
		}

		// Time Complexity: O(log n)
		// insert calls decreaseKey, which has the time complexity of O(log n). See above definition for more details.
		public void insert(HeapNode pointNode)
		{
			int newIndex = count;
			count += 1;
			if (count > maxCountEver) {
				maxCountEver = count;
			}

			if (count > nodes.Count)
			{
				nodes.Add(pointNode);
			}
			else
			{
				nodes[newIndex] = pointNode;
			}

			// bubble up, O(log n)
			decreaseKey(newIndex, pointNode);
		}


		// Time Complexity: O(1)
		private void swap(int index1, int index2)
		{
			// indexMapping[index in originalPoints array] = new index in heap array
			HeapNode temp = nodes[index1];
			nodes[index1] = nodes[index2];
			nodes[index2] = temp;

		}

		public bool isEmpty()
		{
			return count == 0;
		}

	}
}
