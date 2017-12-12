using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;


namespace TSP
{

    public class ProblemAndSolver
    {

        public class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your data structure(s) and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            /// <summary>
            /// constructor
            /// </summary>
            /// <param name="iroute">a (hopefully) valid tour</param>
            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }

            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }

        }

        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        /// <summary>
        /// Default time limit (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Time text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int TIME_LIMIT = 60;        //in seconds

        private const int CITY_ICON_SIZE = 5;


        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        public City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;

        /// <summary>
        /// time limit in milliseconds for state space search
        /// can be used by any solver method to truncate the search and return the BSSF
        /// </summary>
        private int time_limit;
        #endregion

        #region Public members

        /// <summary>
        /// These three constants are used for convenience/clarity in populating and accessing the results array that is passed back to the calling Form
        /// </summary>
        public const int COST = 0;           
        public const int TIME = 1;
        public const int COUNT = 2;
        
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = TIME_LIMIT * 1000;                        // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        public ProblemAndSolver(int seed, int size, int time)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = time*1000;                        // time is entered in the GUI in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// make a new problem with the given size, now including timelimit paremeter that was added to form.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode, int timelimit)
        {
            this._size = size;
            this._mode = mode;
            this.time_limit = timelimit*1000;                                   //convert seconds to milliseconds
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());                    
            else
                return -1D; 
        }

        /// <summary>
        /// This is the entry point for the default solver
        /// which just finds a valid random tour 
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] defaultSolveProblem()
        {
            int i, swap, temp, count=0;
            string[] results = new string[3];
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();
            Stopwatch timer = new Stopwatch();

            timer.Start();

            do
            {
                for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                    perm[i] = i;
                for (i = 0; i < perm.Length; i++)
                {
                    swap = i;
                    while (swap == i)
                        swap = rnd.Next(0, Cities.Length);
                    temp = perm[i];
                    perm[i] = perm[swap];
                    perm[swap] = temp;
                }
                Route.Clear();
                for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                {
                    Route.Add(Cities[perm[i]]);
                }
                bssf = new TSPSolution(Route);
                count++;
            } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"default.txt", true))
            {
                file.WriteLine(Cities.Length.ToString() + "\t" + (((Double)timer.ElapsedMilliseconds) / 1000.0).ToString() + "\t" + costOfBssf().ToString());
            }
            Console.WriteLine("Default finished with {0} cost", results[COST]);
            for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
            {
                Console.WriteLine("Row {0} in path is city {1}", i, Route.IndexOf(Cities[i]));
            }

            return results;
        }


// Time Complexity: O(n^2)
// Space Complexity: O(1)
public double reduceArrayInPlace(ref double[,] costArray) {
    int row = 0, col = 0;
    
    double lowerBound = 0, minValue = 0;
    for (row = 0; row < this._size; row++) {
        minValue = costArray[row, 0];
        for (col = 1; col < this._size; col++) {
            if (costArray[row, col] < minValue) {
                minValue = costArray[row, col];
            }
        }
        if (!double.IsPositiveInfinity(minValue))  {
            // we can reduce
            for (col = 0; col < this._size; col++) {
                if (costArray[row, col] != double.PositiveInfinity) {
                    costArray[row, col] -= minValue;
                }
            }
            lowerBound += minValue;
        }   
    }

    // printCostArray(ref costArray);

    for (col = 0; col < this._size; col++) {
        minValue = costArray[0, col];
        for (row = 1; row < this._size; row++) {
            if (costArray[row, col] < minValue) {
                minValue = costArray[row, col];
            }
        }
        if (minValue < double.PositiveInfinity)  {
            // we can reduce
            for (row = 0; row < this._size; row++) {
                if (costArray[row, col] != double.PositiveInfinity) {
                    costArray[row, col] -= minValue;
                }
            }
            lowerBound += minValue;
        }   
    }
    return lowerBound;
}

// Time Complexity: O(n^2)
// Space Complexity: O(1)
public void fillArrayInPlaceWithCosts(ref double[,] costArray) {
    int row, col = 0;
    for (row=0; row < this._size; row++) {
        for (col=0; col < this._size; col++) {
            if (row == col) {
                costArray[row, col] = double.PositiveInfinity;
            } else {
                costArray[row, col] = Cities[row].costToGetTo(Cities[col]);
            }
        }
    }
}

// helper function that prints a cost array. 
// Time Complexity of O(n^2)
// Space Complexity: O(1)
public void printCostArray(ref double[,] costArray) {
    int row, col = 0;
    double temp = 0;
    for (row=0; row < this._size; row++) {
        for (col=0; col < this._size; col++) {
            temp = costArray[row, col];
            if (temp == double.PositiveInfinity) {
                Console.Write("Inf\t");
            } else {
                Console.Write("{0:N2}\t", temp);
            }
        }
        Console.WriteLine(" ");
    }
}



// Time Complexity: O(n!)
// Space Complexity: O(n^4)

// The time complexity of the Branch and Bound algorithm is O(n!). The algorithm attempts to 
// lower the time required to solve the solutions by pruning the tree of solutions via a 
// lower bound for each partial solution, reducing the number of leaf nodes that are needed to be checked, 
// but in its worst case, the algorithm is still O(n!). This can been seen in a scenario where our algorithm is 
// unable to prune any branches and the BSSF is updated each time our search hits a leaf node. In this case, our algorithm 
// would finally arrive at the optimal cost on the last leaf node, meaning that the algorithm had to deal with every possible 
// solution – in short, the worst case of the B&B algorithm is no different than a brute force approach. Thus, the time complexity 
// is O(n!). Additionally, my implementation of the HeapPriorityQueue allowed all heap operations to be done in O(log n) time. 

public string[] bBSolveProblem()
{
    // run the default solution to get BSSF
    // defaultSolveProblem();            

    // run greedy to speed up :)
    greedySolveProblem();            
    Stopwatch timer = new Stopwatch();

    timer.Start();

    // contains the costs of traveling from any one node to another
    double[,] costArray = new double[this._size, this._size];

    // calculates all the costs and fills our cost array
    fillArrayInPlaceWithCosts(ref costArray);
    double currentLowerBound = reduceArrayInPlace(ref costArray);

    HeapPriorityQueue queue = new HeapPriorityQueue();
    
    ArrayList initialPath;
    int source = 0;

    // add a partial solution starting from each city to our priority queue. 
    // this allows us to check every possible permutation.
    for (source = 0; source < this._size; source++) {
        initialPath = new ArrayList();
        initialPath.Add(source); 
        queue.insert(new HeapNode((double [,]) costArray.Clone(), initialPath, currentLowerBound));
    }

    // temporary variables to prevent repeated variable allocation
    int destination = 0;
    HeapNode currentNode;
    ArrayList currentRoute; 
    double[,] currentCostArray = new double[1, 1];
    double currentBSSFCost = costOfBssf();
    Console.WriteLine("{0} is our currentBSSFCost", currentBSSFCost);
    double costToNextCity = 0;
    double addition = 0;
    bool updatedBSSF = false;
    int updates = 0;
    int prunes = 0;
    int totalInserts = 0;
    bool shouldSkipAfterPath = false;
    TSPSolution tempBssf;

    // Keep working through solutions until our queue is done with partial solutions 
    // or runs out of time.

    // O(n!) time complexity!
    while (timer.Elapsed.TotalMilliseconds < time_limit && !queue.isEmpty()) {
        currentNode = queue.deleteMin();

        // initial prune right after coming off the queue
        if (currentNode.LowerBound >= currentBSSFCost) {
            // less than or equal because once we've got a good path, we don't want anything that's the same
            // Console.WriteLine("pre check pruned node: {0} currentBSSF: {1}", currentNode.LowerBound, currentBSSFCost);
            // PRUNE. 
            prunes++;
            continue;
        }

        // try to add partial solutions visiting all nodes that haven't been visited already
        for (destination = 0; destination < this._size; destination++) {
            costToNextCity = currentNode.costToCity(destination);
            if (double.IsPositiveInfinity(costToNextCity) || currentNode.CurrentPath.Contains(destination)) {
                // skip because we have already visited or should be impossible to visit
            } else {
                totalInserts++;
                shouldSkipAfterPath = false;
                // new state's lower bound starts with the currentNode's lowerBound
                currentLowerBound = currentNode.LowerBound;
                currentLowerBound += costToNextCity;

                if ( currentLowerBound >= currentBSSFCost) {
                    shouldSkipAfterPath = true;
                }

                if (!shouldSkipAfterPath) {
                    currentCostArray = (double [,]) currentNode.CostArray.Clone();
            
                    // no one can come through those paths any more
                    currentCostArray[currentNode.LastCity, destination] = double.PositiveInfinity;

                    // set the cost of returning to a node infinity
                    for (int iter = 0; iter < this._size; iter++) {
                        currentCostArray[currentNode.LastCity, iter] = double.PositiveInfinity;
                        currentCostArray[iter, destination] = double.PositiveInfinity;
                    }

                    // calculate new node lower bound
                    addition = reduceArrayInPlace(ref currentCostArray);   
                    currentLowerBound += addition; 
                }
                            
                // should we prune?
                if (currentLowerBound < currentBSSFCost) {
                    // Console.WriteLine("post check unpruned node: {0} currentBSSF: {1}", currentLowerBound, currentBSSFCost);
                    // no!
                    currentRoute = (ArrayList) currentNode.CurrentPath.Clone();
                    currentRoute.Add(destination);
                    if (currentRoute.Count == this._size) {
                        // we have a complete path and we've beat bssf. update BSSF! 

                        ArrayList newRoute = new ArrayList();
                        foreach (int index in currentRoute) {
                            newRoute.Add(Cities[index]);
                        }
                        // update BSSF with our new route of cities
                        tempBssf = new TSPSolution(newRoute);
                        // calculate the BSSF's actual cost
                        double newBSSFCost = tempBssf.costOfRoute();
                        if (double.IsPositiveInfinity(newBSSFCost)) {
                            Console.WriteLine("WTF this path somehow got a bad edge and the new route is Infinity. This must happen when the next node does not have an edge back to our original node");
                        } else {
                            Console.WriteLine("post check lower bound on node: {0} new bssf: {1} currentBSSF: {2}", currentLowerBound, newBSSFCost, currentBSSFCost);
                            bssf = tempBssf;
                            updatedBSSF = true;
                            currentBSSFCost = newBSSFCost;
                            updates++;
                        }
                    } else {
                        // not yet a complete solution, so just add the partial solution to the state
                        queue.insert(new HeapNode(currentCostArray, currentRoute, currentLowerBound));
                    }
                } else {
                    // Console.WriteLine("post check pruned node: {0} currentBSSF: {1}", currentLowerBound, currentBSSFCost);
                    // PRUNE
                    prunes++;
                }
            }
        }
    }

    timer.Stop();
    string[] results = new string[3];
     
    using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"bnb.txt", true))
            {
                file.WriteLine(Cities.Length.ToString() + "\t" + (((Double)timer.ElapsedMilliseconds)/1000.0).ToString() + "\t" + costOfBssf().ToString());
            }
    results[COST] = costOfBssf().ToString();
    results[TIME] = timer.Elapsed.ToString();

    results[COUNT] = updates.ToString();
    Console.WriteLine("Total Number of States: {0}", totalInserts);
    Console.WriteLine("Most States at Any Point: {0}", queue.maxCountEver);
    Console.WriteLine("Number of Prunes: {0}", prunes);

    return results;
}


        /////////////////////////////////////////////////////////////////////////////////////////////
        // These additional solver methods will be implemented as part of the group project.
        ////////////////////////////////////////////////////////////////////////////////////////////

        /// <summary>
        // finds the greedy tour starting from each city and keeps the best (valid) one
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] greedySolveProblem()
        {

            // workon later!
            string[] results = new string[3];

            double[,] costArray = new double[this._size, this._size];

            // calculates all the costs and fills our cost array
            fillArrayInPlaceWithCosts(ref costArray);

			int i, j, swap, temp, count = 0;
			int[] perm = new int[Cities.Length];
			Route = new ArrayList();
			ArrayList BestRoute = new ArrayList();
            HashSet<int> RouteSet = new HashSet<int>();
            double bestSolutionCost = double.PositiveInfinity;
            double currentBestSolutionCost = double.PositiveInfinity;
			Random rnd = new Random();
			Stopwatch timer = new Stopwatch();

			timer.Start();

            for (i = 0; i < this._size; i++)
            {
                Route.Clear();
                RouteSet.Clear();

                Route.Add(i);
                RouteSet.Add(i);
                currentBestSolutionCost = 0;

                while (RouteSet.Count < this._size)
                {  
                    double minCost = double.PositiveInfinity;
                    int minCity = -1;

                    for (j = 0; j < this._size; j++) 
                    {
                        if (RouteSet.Contains(j)) {
                            continue;
                        } else {
                            if (costArray[ (int) Route[Route.Count - 1], j ] < minCost) {
                                minCost = costArray[ (int) Route[Route.Count - 1], j ];
                                minCity = j;
                            }
                        }
                    }
                    if (minCity < 0) {
                        // no best solution added
                        break;
                    }
                    else {
                        RouteSet.Add(minCity);
                        Route.Add(minCity);
                        currentBestSolutionCost += minCost;
                    }
                }
                if (RouteSet.Count == this._size) {
                    if (currentBestSolutionCost < bestSolutionCost) {
                        bestSolutionCost = currentBestSolutionCost;
                        BestRoute = Route.Clone() as ArrayList;
                        count++;
                    }
                }
            }

            // until a valid route is found
            ArrayList newRoute = new ArrayList();
            foreach (int index in BestRoute) {
                newRoute.Add(Cities[index]);
            }
            // update BSSF with our new route of cities
            bssf = new TSPSolution(newRoute);
            // calculate the BSSF's actual cost
            double newBSSFCost = costOfBssf();

            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
			results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"greedy.txt", true))
            {
                file.WriteLine(Cities.Length.ToString() + "\t" + (((Double)timer.ElapsedMilliseconds) / 1000.0).ToString() + "\t" + costOfBssf().ToString());
            }
            Console.WriteLine("Greedy finished with {0} cost", results[COST]);

            return results;
        }
        public ArrayList Swap(ArrayList path, int i, int k)
        {
            ArrayList newPath = new ArrayList();
            int start = i-1;
            if (start < 0) {
                start = path.Count - 1;
            }
            int end = k+1;
            if (end == path.Count) {
                end = 0;
            }
            // Console.WriteLine("Swapping {0} to {1} with {2} to {3}", i-1, i, k, k+1);

            int x;
            for(x = 0; x < i; x++)
            {
                newPath.Add(path[x]);
            }
            for(x = 0; x <= k-i; x++)
            {
                newPath.Add(path[k-x]);
            }
            for(x = k+1; x < Cities.Length; x++)
            {
                newPath.Add(path[x]);
            }
            ArrayList cities = new ArrayList(Cities);
            // foreach (var c in path) {
            //     Console.Write("{0}->", cities.IndexOf(c));
            // }
            // Console.WriteLine("");
            // foreach (var c in newPath) {
            //     Console.Write("{0}->", cities.IndexOf(c));
            // }
            // Console.WriteLine("");
            return newPath;
        }

        public void recursiveKOpt(ref ArrayList path, int level, ref int[] indexes, ref bool foundSolution) {
            int i = 0;
            if (foundSolution) {
                return;
            }
            if (level == 0) {
                // do array stuff
                Console.Write("(");
                for (i = 0; i < indexes.Length; i++) {
                    Console.Write("{0},", indexes[i]);
                }
                Console.WriteLine(")");

                int temp = 0;
                Double previousBest;
                double newPathCost = 0;
                TSPSolution newPath;

                ArrayList tempPath = path.Clone() as ArrayList;
                for (i = 0; i < indexes.Length - 1; i++) {
                    tempPath = Swap(tempPath, i, i+1);
                }
                newPath = new TSPSolution( tempPath );
                newPathCost = newPath.costOfRoute();

                if(bssf.costOfRoute() > newPathCost)
                {
                    bssf = newPath;
                    foundSolution = true;
                    return;
                }

                for (i = indexes.Length - 1; i >= 0; i--) {
                    indexes[i] += 1;
                    if (indexes[i] == this._size) {
                        indexes[i] = 0;
                    } else {
                        break;
                    }
                }
            } else {
                for (i = 0; i < this._size; i++) {
                    recursiveKOpt(ref path, level-1, ref indexes, ref foundSolution);
                }
            }
        }

        //2-opt 
        public string[] fancySolveProblem()
        {
            Console.WriteLine("--------\nNew K-opt");
            string[] results = new string[3];
            int updates = 0, i, k;
            greedySolveProblem();
            bool betterSolutionFound;
            Double previousBest;
            double newPathCost = 0;
            TSPSolution newPath;
            Stopwatch timer = new Stopwatch();
            timer.Reset();
            timer.Start();

            double[,] costArray = new double[this._size, this._size];

            // calculates all the costs and fills our cost array
            fillArrayInPlaceWithCosts(ref costArray);
            // printCostArray(ref costArray);
            // Console.WriteLine("\n----");
            int level = 2;
            int[] indexes = new int[level];

            for (i=0; i<level; i++) {
                indexes[i] = 0;
            }

            bool foundSolution = false;
            
            recursiveKOpt(ref bssf.Route, level, ref indexes, ref foundSolution);
            
            do
            {
                betterSolutionFound = false;
                recursiveKOpt(ref bssf.Route, level, ref indexes, ref foundSolution);
            } while (betterSolutionFound && timer.Elapsed.TotalMilliseconds < time_limit);

            
            timer.Stop();
            results[COST] = costOfBssf().ToString();
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = updates.ToString();

            return results;
        }

        public string[] fancySolveProblem2()
        {
            Console.WriteLine("--------\nNew K-opt");
            string[] results = new string[3];
            int updates = 0, i, k;
            defaultSolveProblem();
            bool betterSolutionFound;
            Double previousBest;
            double newPathCost = 0;
            TSPSolution newPath;
            Stopwatch timer = new Stopwatch();
            timer.Start();

            double[,] costArray = new double[this._size, this._size];

            // calculates all the costs and fills our cost array
            fillArrayInPlaceWithCosts(ref costArray);
            printCostArray(ref costArray);
            Console.WriteLine("\n----");

            do
            {
                previousBest = costOfBssf();
                betterSolutionFound = false;
                for (i = 0; i < Cities.Length-1 && timer.Elapsed.TotalMilliseconds < time_limit; i++)
                {
                    for(k = i+1; k < Cities.Length-1 && timer.Elapsed.TotalMilliseconds < time_limit; k++)
                    {
                        newPath = new TSPSolution( Swap(bssf.Route, i, k) );
                        newPathCost = newPath.costOfRoute();
                        if(previousBest > newPathCost)
                        {
                            //Console.WriteLine("updating Best after swapping {0} to {1} with {2} to {3}", i-1, i, k, k+1);
                            //Console.WriteLine("Previous best: {0} New Route: {1}", previousBest, newPathCost);
                            bssf = newPath;
                            updates++;
                            betterSolutionFound = true;
                            break;
                        }
                        
                    }
                    if (betterSolutionFound)
                    {
                        break;
                    }
                }
            } while (betterSolutionFound && timer.Elapsed.TotalMilliseconds < time_limit);

            
            timer.Stop();
            results[COST] = costOfBssf().ToString();
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = updates.ToString();
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"special.txt", true))
            {
                file.WriteLine(Cities.Length.ToString() + "\t" + (((Double)timer.ElapsedMilliseconds) / 1000.0).ToString() + "\t" + costOfBssf().ToString());
            }
            return results;
        }
        #endregion
    }

}
