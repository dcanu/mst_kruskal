/**
 * File: mst.java
 * Implements Kruskal algorithm, using subroutines of Quick Sort and
 * Directed Forest for disjoint set.
 * <p>
 * COMP6466 Assignment 2
 * Student ID: u4610248 Debashish Chakraborty
 * <p>
 * The file can be compiled in a standard method from command line,
 * since it does not take any input from command line
 * javac mst.java
 * java mst.class
 * <p>
 * The variations required for execution with different random weights are:
 * To change between random weights between [0,1] (default when submitted) and [0,0.5]
 * Step 1: Change the rndRange range (commented) in the generateGraph method
 * Step 2: Change the print out (commented) in the main method
 */


import java.util.*;

public class mst {
    // represent the mst as a list of edges
    ArrayList<Edge> mstree = new ArrayList<>();

    // weight of the MST and a getter method to retrieve it
    double weight;

    double getWeight() {
        return weight;
    }

    /*
     *  The main method prints out average weighted sum for 15 MSTs and the average
     *  running time of Kruskal's algorithm for finding the MSTs
     */
    public static void main(String[] args) {

        if (args.length != 0)
            throw new InputMismatchException("No input " +
                    "required for the MST");

        // given list of vertex size
        int[] numVlist = {10, 100, 200, 500, 1000};

        for (int i = 0; i < numVlist.length; i++) {

            int nVertices = numVlist[i];

            double[] avgWeight = new double[numVlist.length];
            long[] avgTime = new long[numVlist.length];

            // creating 15 graphs
            for (int j = 0; j < 15; j++) {

                double[][] graph = new double[nVertices][nVertices];
                mst.generateGraph(graph);

                long startTime = System.nanoTime();

                mst mst = new mst(graph);
                avgWeight[i] += mst.getWeight() / 15;
                long endTime = System.nanoTime();

                avgTime[i] += ((endTime - startTime) / 15);

            }

//            // print one line of output for each n
            System.out.printf("Vertices: %d | avg.L(n): %s units | avg.T(n): %dns%n",
                    nVertices, avgWeight[i], avgTime[i]);

            // printout when the weight is random number [0,0.5]
//            System.out.printf("Vertices: %d | avg.L'(n): %s units | avg.T'(n): %dns%n",
//                    nVertices, avgWeight[i], avgTime[i]);


        }

    }

    /*
    *   Create the MST taking the graph as input using Quicksort to sort the edges
    *   according to its corresponding weights and Kruskal algorithm using disjoint
    *   forest structure to connect vertices in order to avoid cycle
    *
     */

    public mst(double[][] G) {

        int graphSize = G.length;
        // number of edges of a complete graph
        int numEdges = graphSize * (graphSize - 1) / 2;
        Comparable[] edges = new Comparable[numEdges];
        // counter for edges
        int count = 0;
        // create a list of edges from the graph
        for (int i = 0; i < graphSize; i++) {
            for (int j = i; j < graphSize; j++) {
                if (i == j)
                    continue;
                else {
                    Edge edge = new Edge(i, j, G[i][j]);
                    edges[count++] = edge;
                }
            }
        }

        QuickSort.sort(edges);

        DisjointSet ds = new DisjointSet(graphSize);
        // count number of edges remaining
        int count_rem = edges.length;

        while ((count_rem > 0) && mstree.size() < graphSize - 1) {
            Edge edge = (Edge) edges[count - count_rem];
            count_rem--;
            int u = edge.getU();
            int v = edge.getV(u);
            // check if they already have the same root
            // if vertices do not have the same root, using disjoint set union
            // and add the edge to the mst
            if (!ds.isConnected(u, v)) {

                ds.union(u, v);
                mstree.add(edge);
                weight += edge.getWeight();

            }
        }
    }

    // Store the graph as a symmetric matrix containing the graph edges
    // This includes both instances where graph needs to
    public static double[][] generateGraph(double[][] graph) {

        int size = graph.length;

        for (int i = 0; i < size; i++) {
            for (int j = i; j < size; j++) {
                if (i == j) {
                    graph[i][j] = 0;
                } else {
                    // for weights between [0,1]
                    graph[i][j] = rndRange(0, 1);
                    // for weights between [0,0.5]
//                    graph[i][j] = rndRange(0.0, 0.5);
                    graph[j][i] = graph[i][j];
                }
            }
        }

        return graph;

    }

    // print out random double for a designated range
    public static double rndRange(double min, double max) {
        double range = Math.abs(max - min);
        return (Math.random() * range) + (min <= max ? min : max);
    }

    static class QuickSort {
        static Random random = new Random();

        QuickSort() {
        }

        /**
         * Sort the array elements in a non-decreasing order
         */
        static void sort(Comparable[] a) {
            shuffle(a);
            sort(a, 0, a.length - 1);
        }

        /*
         * sort the subarray from a[left] to a[right]
          */
        static void sort(Comparable[] a, int left, int right) {
            if (right <= left)
                return;
            int j = partition(a, left, right);
            sort(a, left, j - 1);
            sort(a, j + 1, right);
        }

        /*
        *  Arbitrarily choose a[left] to be the partitioning item - the one that will go to the
        *  final position. Then scan from the left end of the array until we find an entry that
        *  is greater than (or equal to) the partitioning item, and scan the right end of the array
        *  until we find an entry less than (or equal to) the partitioning item.
        *   partition the subarray a[left..right] such that
        *    a[left..j-1] <= a[j] <= a[j+1..right]
        *    and return the index j.
         */
        static int partition(Comparable[] a, int first, int last) {
            int i = first;
            int j = last + 1;
            Comparable v = a[first];

            while (true) {
                while (a[++i].compareTo(v) < 0)
                    if (i == last)
                        break;
                while (v.compareTo(a[--j]) < 0)
                    if (j == first)
                        break;
                if (i >= j)
                    break;
                swap(a, i, j);
            }
            swap(a, first, j);
            return j;
        }


        // swap a[i] and a[j]
        static void swap(Object[] a, int i, int j) {
            Object temp = a[i];
            a[i] = a[j];
            a[j] = temp;
        }

        /**
         * Rearrange the elements of an array in random order
         * between i and N-1
         */
        static void shuffle(Object[] a) {

            int N = a.length;
            for (int i = 0; i < N; i++) {
                int r = i + randomize(N - i);
                Object temp = a[i];
                a[i] = a[r];
                a[r] = temp;
            }
        }

        /**
         * Returns an integer uniformly between 0 (inclusive) and N (exclusive)
         */
        static int randomize(int N) {
            return random.nextInt(N);
        }
    }

    /*
    *   In the disjoint set forest each member points only to its parent. The
    *   root of each tree contains the representative and is its own parent.
     */

    static class DisjointSet {

        int[] id, rank;

        /*
        *   make_set operation -
        *   creates a new set whose only member (and thus representative)
        *   is itself
         */
        public DisjointSet(int size) {

            id = new int[size];
            rank = new int[size];

            for (int i = 0; i < size; i++) {
                id[i] = i;
                rank[i] = 0;
            }
        }

        /*
        *   find_set returns a pointer to the representative
        *   of the (unique) set containing the node
         */

        public int find_set(int n) {

            while (n != id[n])
                n = find_set(id[n]);

            return n;
        }

        /*
        *   Link uses is used in union to assign parents and link
        *   trees into forest according to rank hueristic
         */

        public void link(int x, int y) {

            if (rank[x] > rank[y]) {
                id[y] = x;
            } else {
                id[x] = y;
            }

            if (rank[x] == rank[y]) {
                rank[y]++;
            }

        }

        /*
        *   union uses path compression to join roots of disjoint set tree using link
        *   for path compression, thus creating an efficient disjoint set forest
        *
         */

        public void union(int x, int y) {

            int id_x = find_set(x);
            int id_y = find_set(y);


            if (id_x == id_y) {
                return;
            }

            link(id_x, id_y);

        }


        /*
        *   Check if the roots are the same for two vertices
         */
        public boolean isConnected(int x, int y) {
            return find_set(x) == find_set(y);
        }


    }

    /*
    *   Edge is defined by its vertices and weight respectively
    *   We are also able to distinguish between two edges of an
    *   undirected graph, which necessarily represent the same edge
    *   such as (1,3) and (3,1).
    *   Here, such edges are only counted once using getU and getV methods
     */

    static class Edge implements Comparable<Edge> {

        int u, v;
        double w;

        public Edge(int u, int v, double w) {
            this.u = u;
            this.v = v;
            this.w = w;
        }

        public double getWeight() {
            return w;
        }

        public int getU() {
            return u;
        }

        /*
        During the implementation of disjoint set, it is crucial not to consider the same edge twice
        since the graph is undirected. Edge (1,2) and (2,1) should be counted as the same edge.
         */

        public int getV(int w) {
            if (w == u)
                return v;
            else
                return u;
        }

        @Override
        public int compareTo(Edge e) {
            if (!(e instanceof Edge)) {
                throw new ClassCastException("Not an instance of Edge");
            }

            if (this.getWeight() < e.getWeight()) {
                return -1;
            }

            if (this.getWeight() > e.getWeight()) {
                return 1;
            }

            if (this == e) {
                return 0;
            }

            return 0;
        }
    }


}
