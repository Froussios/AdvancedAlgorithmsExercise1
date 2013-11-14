//
//  greedy.java
//  AdvancedProject
//
//  Created by Sicco on 9/17/09.
//  Copyright 2009 __MyCompanyName__. All rights reserved.
//

// for file reading and writing
import java.io.*;
import java.util.Scanner;

import java.util.ArrayList;
// for the best-first (and breadth-first) search
import java.util.PriorityQueue;
import java.util.Comparator;

// can be useful to create an efficient representation of a schedule
import java.util.BitSet;
// can be usefull for memoization
import java.util.HashMap;

// contains usefull functions
import java.lang.Math;

import java.util.*; 

import javax.swing.event.ListSelectionEvent;

// The schedule class, represents a (partial) sequence of jobs
class schedule implements Comparable {
	
	// A linked-list is a relatively efficient representation of a schedule
	// Feel free to modify it if you feel there exists a better one
	// The main advantage is that in a search-tree there is a lot of overlap
	// between schedules, this implementation stores this overlap only once
	public schedule previous;
	public int scheduled_job;
	
	// tardiness can be calculated instead of memorised
	// however, we need to calculate it a lot, so we memorise it
	// if memory is an issue, however, try calculating it
	public int tardiness;
	
	public schedule(){
		scheduled_job = -1;
		previous = null;
		tardiness = 0;
	}
	
	// add an additional job to the schedule
	public schedule(schedule s, int job){
		previous = s;
		scheduled_job = job;
		tardiness = Math.max(0, get_total_time() - algorithms.jobs[scheduled_job][1]);
		if(previous != null)
			tardiness += previous.tardiness;
	}
	
	// used by the best-first search
	// currently, schedules are traversed in smallest total tardiness order
	public int compareTo(Object o){
		return (get_tardiness()) - (((schedule)o).get_tardiness());
		
		// replace with the following to get a depth-first search
		// return get_depth() - ((schedule)o).get_depth();
	}
	
	public int get_depth(){
		if(previous != null)
			return previous.get_depth() + 1;
		return 1;
	}
	
	public int get_total_time(){
		if(previous != null)
			return previous.get_total_time() + algorithms.jobs[scheduled_job][0];
		return algorithms.jobs[scheduled_job][0];
	}
	
	public int get_tardiness(){
		return tardiness;
	}
	
	public boolean contains(int job){
		return (scheduled_job == job) || (previous != null && previous.contains(job));
	}
}

class greedy {
	// returns the earliest deadline first schedule
	// sorting is a little quicker, but then it is more tricky
	// to use this as a subroutine for a search method
	public static schedule greedy(){
		int due = -1;
		int job_to_schedule = -1;
		for(int i = 0; i < algorithms.num_jobs; ++i){
			if(due == -1 || due > algorithms.jobs[i][1]){
				due = algorithms.jobs[i][1];
				job_to_schedule = i;
			}
		}
		return greedy(new schedule(null, job_to_schedule));
	}
	
	// adds the next earliest deadline first job to the schedule
	public static schedule greedy(schedule s){
		if(s.get_depth() >= algorithms.num_jobs) return s;
		
		int due = -1;
		int job_to_schedule = -1;
		for(int i = 0; i < algorithms.num_jobs; ++i){
			if(s.contains(i) == false && (due == -1 || due > algorithms.jobs[i][1])){
				due = algorithms.jobs[i][1];
				job_to_schedule = i;
			}
		}
		
		s = new schedule(s, job_to_schedule);
		return greedy(s);
	}
}

class best_first_search {
	// returns the best-first (or breadth-first) search schedule
	// It uses a PriorityQueue to store schedules, in every iteration
	// it gets the next best schedule, tries to append all possible jobs
	// and stores the resulting schedules on the queue
	public static schedule search(){
		PriorityQueue<schedule> Q = new PriorityQueue<schedule>();
		
		for(int i = 0; i < algorithms.num_jobs; ++i){
			Q.offer(new schedule(null, i));
		}
		
		schedule best_schedule = null;
		while(Q.peek() != null){
			schedule s = Q.poll();
			
			// can be useful for debugging
			// System.err.println(Q.size());

			if(s.get_depth() >= algorithms.num_jobs){
				if(best_schedule == null || best_schedule.get_tardiness() > s.get_tardiness()){
					best_schedule = s;
				}
				continue;
			}
			
			if(best_schedule != null && best_schedule.get_tardiness() < s.get_tardiness()){
				continue;
			}
			
			for(int i = 0; i < algorithms.num_jobs; ++i){
				if(s.contains(i) == false){
					Q.offer(new schedule(s, i));
				}
			}
		}
		return best_schedule;
	}
}

class Dynamic
{	
	private static int[][][][] memory = null;
	private static void initializeMemory()
	{
		int n = algorithms.jobs.length;
		int P = 0;
		for (int i = 0; i < algorithms.jobs.length; i++)
			P += Job.length(i);
		
		memory = new int[n][n][n][P];
		// Initialise all the elements to a non-value
		for (int i1=0 ; i1<memory.length ; i1++ )
			for (int i2 = 0; i2<memory[i1].length ; i2++)
				for (int i3 = 0; i3 < memory[i1][i2].length; i3++)
					for (int i4 = 0; i4 < memory[i1][i2][i3].length; i4++)
						memory[i1][i2][i3][i4] = -3;
	}
	
	public static int searchRec()
	{
		Subset problem = new Subset();
		
		initializeMemory();
		
		ArrayList<Integer> sequence = new ArrayList<Integer>();
		int tardiness = recursive(problem, 0, 0, sequence);
		
//		if ( sequence.size() != algorithms.num_jobs )
//		{
//			System.out.println("(Missing jobs in solution) ");
//		}
//		int real_tard = 0;
//		int time = 0;
//		for (Integer i : sequence)
//		{
//			System.out.print("" + i + "-");
//			time += Job.length(i);
//			real_tard += Math.max(0, time - Job.deadline(i));
//		} 
//		System.out.print(">" + real_tard);
//		System.out.println();
//		System.out.println("Hits: " + hits + " Misses: " + misses);
		
		return tardiness;
	}
	
	public static int search()
	{
		initializeMemory();
		
		Subset problem = new Subset();
		
		int tardiness = run(problem);
		
		return tardiness;
	}
	
	static int run(Subset problem)
	{
		for ( int i=0 ; i<memory.length ; i++ )
		{
			for ( int j=i ; j<memory[i].length ; j++ )
			{
				for ( int k=0 ; k<=j ; k++ ) //TODO: k<=j ok?
				{
					for ( int t=0 ; t<memory[i][j][k].length ; t++ ) //TODO: can this be limited to P - sum(length(a): a in S(i,j,k))
					{
						int minTard = Integer.MAX_VALUE;
						
						ExperimentalSubset set = new ExperimentalSubset(i,j,k);
						ArrayList<Integer> jobs = set.allElements();
						
						if ( jobs.size() == 0 )
						{
							// Set has no items
							minTard = 0;
						}
						else if ( jobs.size() == 1 )
						{
							// Set has exactly one item
							int job = jobs.get(0);
							minTard = Math.max(0, t + Job.length(job) - Job.deadline(job));
						}
						else
						{
							int pivot = set.maxLength();
							ArrayList<Integer> subset1 = new ArrayList<Integer>(jobs.subList(0, jobs.indexOf(pivot)));
							ArrayList<Integer> subset2 = new ArrayList<Integer>(jobs.subList(jobs.indexOf(pivot)+1, jobs.size()-1));
							
							while (true)
							{
								ExperimentalSubset ss1 = new ExperimentalSubset(subset1);
								ExperimentalSubset ss2 = new ExperimentalSubset(subset2);
								
								int pivotCompletion = t + ss1.totalLength();
								int tard1 = memory[ss1.i][ss1.j][ss1.k][t];
								int tard2 = memory[ss2.i][ss2.j][ss2.k][t+pivotCompletion];
								int tardPivot = Math.max(0, pivotCompletion - Job.deadline(pivot));
								
								int totalTardiness = tard1 + tardPivot + tard2;
								
								if ( totalTardiness < minTard )
									minTard = totalTardiness;
							}
						}
						
						memory[i][j][k][t] = minTard;
					}
				}
			}
		}
		return -1;
	}
	
	static int hits = 0;
	static int misses = 0;
	// Pseudo-recursive implementation of dynamic solution
	static int recursive(Subset subset, int t, int depth, ArrayList<Integer> rvSeq)
	{
		// Count items
		int n= subset.count();
		
		// Subset is empty
		if (n == 0)
			return 0;
		
		// Check if solution already found
		if ( subset.k >- 1 )
		{
			if ( memory[subset.i][subset.j][subset.k][t] >=0 )
			{
				hits++;
				return memory[subset.i][subset.j][subset.k][t];
			}
		}
		misses++;
		
		// Find k : k belongs to subset and has maximum length
		int k = subset.maxLength();
		
		int best = Integer.MAX_VALUE;
		//ArrayList<Integer> best_sequence = null;
		
		for ( int split=k ; split <= subset.j ; split++ )
		{
			// Calculate subsets
			if ( subset.k!=-1 )
				if ( Job.length(subset.k) <= Job.length(k) )
					System.out.print("");
			
			// TODO: use subset.k for subsets contain items of length k
			Subset subset1 = new Subset(subset.i,split,k, subset);
			Subset subset2 = new Subset(split+1,subset.j,k, subset);
			
			ArrayList<Integer> seqPre = null;// new ArrayList<Integer>();
			ArrayList<Integer> seqPost = null; //new ArrayList<Integer>();
			
			int completion_k = t + subset1.totalLength() + Job.length(k);
			
			// Calculate tardiness for subsets, k and total tardiness
			int tard1 = recursive(subset1, t, depth+1, seqPre);
			int tard2 = recursive(subset2, completion_k, depth+1, seqPost);
			int tardk = Math.max(0, completion_k - Job.deadline(k));
			
			int tard = tard1 + tardk + tard2;
			
			// If new optimal solution is found
			if ( tard < best )
			{
				best = tard;
				//best_sequence = seqPre;
				//best_sequence.add(k);
				//best_sequence.addAll(seqPost);
			}
		}
		
		if ( best == Integer.MAX_VALUE )
			System.err.println("Failed to find solution: S" + subset + " k:" + k);
		
		//rvSeq.addAll(best_sequence);
		
		if ( subset.k > -1 )
			memory[subset.i][subset.j][subset.k][t] = best;
		
		return best;
	}
}

class Recursive
{		
	public static int search()
	{
		Problem problem = new Problem(algorithms.num_jobs);
		for ( int j=0 ; j<algorithms.num_jobs ; j++ )
			problem.add(new Job(j));
		
		
		ArrayList<Integer> sequence = new ArrayList<Integer>();
		int result = recursive(problem, 0, sequence);
		
//		int real_tard = 0;
//		int time = 0;
//		for (Integer i : sequence)
//		{
//			System.out.print("" + i + "-");
//			time += Job.length(i);
//			real_tard += Math.max(0, time - Job.deadline(i));
//		} 
//		System.out.print(">" + real_tard);
//		System.out.println();
		
		return result;
	}
	
	static int recursive(Problem problem, int start, ArrayList<Integer> rvSeq)
	{
		if (problem.size() == 0)
			return 0;
		
		int n = problem.size();
		int k = problem.maxLength();
		
		int best = Integer.MAX_VALUE;
		//ArrayList<Integer> seqBest = null;
		
		for ( int d=0 ; d <= n-k-1 ; d++ )
		{
			Problem subproblem1 = new Problem(k+d/*+1*/);
			for ( int i=0 ; i<k+d+1 ; i++ )
				if ( i!=k )
					subproblem1.add(problem.get(i));
			
			Problem subproblem2 = new Problem(n - (k+d));
			for ( int i=k+d+1 ; i<n ; i++ )
				subproblem2.add(problem.get(i));
			
			ArrayList<Integer> seqPre = null;//new ArrayList<Integer>();
			ArrayList<Integer> seqPost = null;//new ArrayList<Integer>();
			
			int tard1 = recursive(subproblem1, start, seqPre);
			int finish_k = start + subproblem1.totalLength() + problem.get(k).length;
			int tard2 = recursive(subproblem2, finish_k, seqPost);
			
			int tard_k = tard1 + tard2 + Math.max(finish_k-problem.get(k).deadline, 0);
			
			if (tard_k < best)
			{
				best = tard_k;
				//seqBest = seqPre;
				//seqBest.add(problem.get(k).job);
				//seqBest.addAll(seqPost);
			}
		}
		
		//System.out.println("\tJob " + problem.get(k) + " start: " + best_preceeding.totalLength() + " tardiness " + Math.max(0, start - problem.get(k).deadline));
		
		//rvSeq.addAll(seqBest);
		
		return best;
	}
}

class BruteForce
{	
	public static int search()
	{
		Problem problem = Problem.LoadProblem();
		
		// Try all permutations
		ArrayList<ArrayList<Job>> permutations = new ArrayList<ArrayList<Job>>();
		permute(new ArrayList<Job>(), problem, permutations);
		
		int minTard = Integer.MAX_VALUE;
		ArrayList<Job> best_sequence = null;
		for (ArrayList<Job> sequence : permutations)
		{
			int time = 0;
			int totalTard = 0;
			for (Job job : sequence)
			{
				time += job.length;
				totalTard += Math.max(0, time - job.deadline);
			}
			if ( totalTard < minTard)
			{
				minTard = totalTard;
				best_sequence = sequence;
			}
		}
		return minTard;
	}
	
	static void permute(ArrayList<Job> arranged, ArrayList<Job> remaining, ArrayList<ArrayList<Job>> permutations)
	{
		if (remaining.size() == 0)
			permutations.add(arranged);
		
		for (Job job : remaining)
		{
			ArrayList<Job> cloneRemaining = (ArrayList<Job>)remaining.clone();
			ArrayList<Job> cloneArranged = (ArrayList<Job>)arranged.clone();
			
			cloneArranged.add(job);
			cloneRemaining.remove(job);
			
			permute(cloneArranged, cloneRemaining, permutations);
		}
	}
}

class ExperimentalSubset
{	
	int i,j,k;
	
	public ExperimentalSubset(int i, int j, int k)
	{
		this.i = i;
		this.j = j;
		this.k = k;		
	}
	
	// This receives a collection of elements and reverse-engineers
	// and reverse-engineers an appropriate definition of the form S(i,j,k)
	public ExperimentalSubset(List<Integer> sortedElements)
	{
		// Sort by job id (already in non-decreasing deadline order)
		//Collections.sort(sortedElements); // Given?
		this.i = sortedElements.get(0);
		this.j = sortedElements.get(sortedElements.size()-1);
		
		int maxElementLength = Integer.MIN_VALUE;
		for ( Integer i : sortedElements )
			if ( maxElementLength < Job.length(i) )
				maxElementLength = Job.length(i);
		
		// Scan the entire problem to find a k that is minimally bigger than maxElementLength
		this.k = -1;
		int kVal = Integer.MAX_VALUE;
		for ( int i=0 ; i<algorithms.num_jobs ; i++ )
			if ( Job.length(i) - maxElementLength > 0 && Job.length(i) < kVal )
			{
				this.k = i;
				kVal = Job.length(i);
			}
		if ( kVal == Integer.MAX_VALUE )
			// No k was found with length greater than every element in the element list
			// The element list contains one of the elements with the maximum length possible
			kVal += 0; //TODO: handle it
	}
	
	public boolean contains(int job)
	{
		if ( k==-1 )
			return job >= i && job <= j ;
		else
			return job >= i 
				&& job <= j 
				&& Job.length(job) < Job.length(k);
	}
	
	public int count()
	{
		int n=0;
		for ( int i=this.i ; i<=this.j ; i++ )
			if ( this.contains(i) )
				n++;
		return n;
	}
	
	public int maxLength()
	{
		int k = -22;
		for ( int i=this.i ; i<=this.j ; i++ )
			if ( this.contains(i) )
				if ( k==-22 )
					k = i;
				else
					if ( Job.length(k) < Job.length(i))
						k = i;
		return k;
	}
	
	public int totalLength()
	{
		int total = 0;
		for ( int i=this.i ; i<=this.j ; i++ )
			if ( this.contains(i) )
				total += Job.length(i);
		return total;
	}
	
	public ArrayList<Integer> allElements()
	{
		ArrayList<Integer> rv = new ArrayList<Integer>();
		for ( int job = i ; job <= j ; job++ )
			if (  contains(job) )
				rv.add(job);
		return rv;
	}
	
	public String toString()
	{
		return "S("+i+","+j+","+k+")";
	}
}

@SuppressWarnings("serial")
class Subset extends ArrayList<Integer>
{
	public int i,j,k;
	
	// Loads the full problem
	public Subset()
	{
		super(algorithms.num_jobs);
		this.i = 0;
		this.j = algorithms.num_jobs-1;
		this.k = -1;
		for (int i = 0; i < algorithms.num_jobs; i++) 
			this.add(i);
		
	}
	public Subset(int start/*inclusive*/, int end/*inclusive*/, int pivot, Subset parentSet)
	{
		super(parentSet.size());
		this.i = start;
		this.j = end;
		this.k = pivot;
		for ( Integer i : parentSet )
			if ( i >= start && i <= end )
				if ( i != pivot )
					this.add(i);
	}
	
//	public boolean contains(int job)
//	{
//		if ( k==-1 )
//			return job >= i && job <= j ;
//		else
//			return job >= i 
//				&& job <= j 
//				&& Job.length(job) < Job.length(k);
//	}
	
	public int count()
	{
//		int n=0;
//		for ( int i=this.i ; i<=this.j ; i++ )
//			if ( this.contains(i) )
//				n++;
//		return n;
		return this.size();
	}
	
	public int maxLength()
	{
//		int k = -22;
//		for ( int i=this.i ; i<=this.j ; i++ )
//			if ( this.contains(i) )
//				if ( k==-22 )
//					k = i;
//				else
//					if ( Job.length(k) < Job.length(i))
//						k = i;
//		return k;
		int k = -22;
		for ( Integer i : this )
			if ( k < 0 )
				k = i;
			else if ( Job.length(k) < Job.length(i) )
				k = i;
		return k;
	}
	
	public int totalLength()
	{
		int total = 0;
//		for ( int i=this.i ; i<=this.j ; i++ )
//			if ( this.contains(i) )
//				total += Job.length(i);
		for ( Integer i : this )
			total += Job.length(i);
		return total;
	}
	
	public String toString()
	{
		return "S(" + i + "," + j + "," + k + ")";
	}
	public boolean equals(Subset s2)
	{
		return this.i == s2.i
				&& this.j == s2.j
				&& this.k == s2.k;
	}
}

class Job
{
	static int length(int job) { return algorithms.jobs[job][0]; }
	static int deadline(int job) { return algorithms.jobs[job][1]; }
	
	public final int job;
	public final int length;
	public final int deadline;
	public Job(int j)
	{
		this.job = j;
		this.length = length(j);
		this.deadline = deadline(j);
	}
	
	public String toString()
	{
		return "" + job + "(" +length+","+deadline+")";
	}
}

@SuppressWarnings("serial")
class Problem extends ArrayList<Job>
{
	public Problem(int size)
	{
		super(size);
	}
	public Problem(int start, int end, int k)
	{
		Job jobk = new Job(k);
		for ( int i=start ; i<=end ; i++ )
		{
			Job job = new Job(i);
			if ( job.length < jobk.length )
				this.add(job);
		}
	}
	
	public static Problem LoadProblem()
	{
		Problem p = new Problem(algorithms.num_jobs);
		for ( int j=0 ; j<algorithms.num_jobs ; j++ )
			p.add(new Job(j));
		return p;
	}
	
	// Returns the index of the job with the maximum processing time
	public int maxLength()
	{
		// Guarantee that Subproblem is not empty?
		int max = 0;
		for ( int j = 0 ; j < this.size() ; j++ )
			if (get(j).length > get(max).length)
				max = j;
		return max;
	}
	
	// Returns the total length of all the jobs in the problem
	public int totalLength()
	{
		int sum = 0;
		for (Job job : this)
			sum += job.length;
		return sum;
	}
}

class algorithms {
	static int num_jobs;
	// size = [num_jobs][2], for every job [0] is the length, [1] is the due time
	static int jobs[][];
	
	// reading a minimum tardiness scheduling problem from a file
	public static void read_problem(String text_file){
		Scanner s = null;
		try {
			s = new Scanner(new BufferedReader(new FileReader(text_file)));
			if(s.hasNextInt()){
				num_jobs = s.nextInt();
				jobs = new int[num_jobs][2];
				int job = 0;
			
				while (s.hasNextInt() && job < num_jobs) {
					jobs[job][0] = s.nextInt();
					jobs[job][1] = s.nextInt();
					job++;
				}
			}
			s.close();
		} catch(Exception e) {
			System.err.println(e);
		}
	}

	public static void main(String args[])
	{
		read_problem(args[0]);
		
		// Sort problem in non-decreasing order of deadline
		Arrays.sort(algorithms.jobs, new Comparator<int[]>() {
				public int compare(int[] f1, int[] f2)
					{
		            	return f1[1] - f2[1];
					}
			});
		
		
		{
			// DEBUGGING: Print problem
			System.out.println("Problem: id(length,due)");
			for (int i=0 ; i < algorithms.num_jobs ; i++ )
				System.out.print(new Job(i) + " , ");
			System.out.println();
		}
		
		interactiveDebugging();
		
		runAlgorithms();
	}

	private static void interactiveDebugging() {
		String strin = "";
		Scanner in = new Scanner(new InputStreamReader(System.in));
		while ( !in.equals("done") )
		{
			strin = in.nextLine();
			String[] elements = strin.split(",");
			ArrayList<Integer> ints = new ArrayList<Integer>();
			for (String str : elements)
				ints.add(Integer.parseInt(str));
			ExperimentalSubset es = new ExperimentalSubset(ints);
			System.out.println(es);
		}
	}
	
	// reads a problem, and outputs the result of both greedy and best-first
    public static void runAlgorithms ()
    {		
		// Run brute force version of dynamic solution
		// Only for the smallest problems
		//System.out.print("Brute force: ");
		//System.out.println(BruteForce.search());
		
		// Run dynamic programming solution (recursive method for debugging)
		System.out.print("Dynamic (Recursive): ");
		System.out.println(Dynamic.searchRec());
		
		// Run dynamic programming solution
		System.out.print("Dynamic: ");
		System.out.println(Dynamic.search());
		
		// Run recursive version of dynamic solution
		System.out.print("Recursive: ");
		System.out.println(Recursive.search());
		
		// Run greedy solution
		System.out.print("Greedy: ");
		schedule s = greedy.greedy();
		System.out.println(s.get_tardiness());
		
		// Run best first solution
		boolean inTheMoodForCrashing = true;
		if (inTheMoodForCrashing)
		{
			try {
				System.out.print("Best first: ");
				s = best_first_search.search();
				System.out.println(s.get_tardiness());
			} catch(Throwable e) {
				// catches out of memory errors
				//e.printStackTrace();
				System.out.println("crashed");
			}
		}
	}
}
