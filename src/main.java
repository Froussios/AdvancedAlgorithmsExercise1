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
	// the memoy used by the dynamic pogramming implementation
	private static int[][][][] memory = null;
	private static void initializeMemory()
	{
		int n = algorithms.jobs.length;
		int P = 0;
		for (int i = 0; i < algorithms.jobs.length; i++)
			P += Job.length(i);
		
		memory = new int[n][n][n+1][P+1];
		// Initialise all the elements to a non-value
		for (int i1=0 ; i1<memory.length ; i1++ )
			for (int i2 = 0; i2<memory[i1].length ; i2++)
				for (int i3 = 0; i3 < memory[i1][i2].length; i3++)
					for (int i4 = 0; i4 < memory[i1][i2][i3].length; i4++)
						memory[i1][i2][i3][i4] = Integer.MAX_VALUE;
	}
	
	// Runs the iterative implementation of the dynamic programming solution
	public static int search()
	{
		initializeMemory();
		
		int tardiness = iterative();
		
		return tardiness;
	}
	
	// Creates a sorted list of all the jobs in the specified range.
	// The items are in non-decreasing order of length
	public static ArrayList<Integer> sortByLength(int rangeStart, int rangeEnd)
	{
		ArrayList<Integer> pivotsInOrderOfLength = new ArrayList<Integer>();
		for ( int k=rangeStart ; k<=rangeEnd ; k++ ) 
			pivotsInOrderOfLength.add(k);
		
		Collections.sort(pivotsInOrderOfLength, new Comparator<Integer>() {
			public int compare(Integer j1, Integer j2)
			{
				return Job.length(j1) - Job.length(j2);
			}
		});
		
		return pivotsInOrderOfLength;
	}
	
	static int iterative()
	{		
		int n = algorithms.num_jobs;
		ArrayList<Integer> jobsByLength = sortByLength(0, n-1);
		jobsByLength.add(n);
		
		for ( int frame = 0 ; frame < n ; frame++ )
		{
			for ( int i=0 ; i < n-frame ; i++ )
			{
				int j = i + frame;				
				for ( Integer setk : jobsByLength )
				{
					if ( setk > j  && setk!=n )
						continue;
					
					Subset set = (setk == n) ? new Subset(i,j,-1) : new Subset(i,j,setk);
					
					for ( int t=0 ; t<memory[i][j][setk].length ; t++ )
					{
						int minTard = Integer.MAX_VALUE;
						int setSize = set.count();
						
						if ( setSize == 0 )
						{	// Set has no items
							minTard = 0;
						}
						else if ( setSize == 1 )
						{	// Set has exactly one item
							int job = set.first();
							minTard = Math.max(0, t + Job.length(job) - Job.deadline(job));
						}
						else
						{
							int pivot = set.maxLength();
							
							for ( int split=pivot ; split <= set.j ; split++ )
							{
								// Calculate subsets
								Subset s1 = new Subset(set.i, split, pivot);
								Subset s2 = new Subset(split+1,set.j, pivot);
								
								int pivotCompletion = t + s1.totalLength() + Job.length(pivot);
								if ( pivotCompletion >= memory[0][0][0].length )
									break;
								
								// Calculate tardiness for subsets, k and total tardiness
								int tard1 = memory[s1.i][s1.j][s1.k][t];
								int tard2 = (s2.i <= s2.j) ? memory[s2.i][s2.j][s2.k][pivotCompletion] : 0;
								int tardk = Math.max(0, pivotCompletion - Job.deadline(pivot));								
								
								int tard = tard1 + tardk + tard2;
								if ( tard1 == Integer.MAX_VALUE || tard2 == Integer.MAX_VALUE )
									tard = Integer.MAX_VALUE;
								
								if ( tard < minTard ) // Improved solution found
									minTard = tard;
							}
						}
						
						memory[i][j][setk][t] = minTard;
					}
				}
			}
		}

		return memory[0][algorithms.num_jobs-1][algorithms.num_jobs][0];
	}
}


// Implements the subset S(i,j,k), as defined in section 3 of the Lawler (1977) paper.
// The definition is extended to include jobs that have processing time equal to k and
// are positioned before k in a non-decreasing deadline ordering of the jobs
class Subset
{	
	public int i,j,k;
	
	public Subset(int i, int j, int k)
	{
		this.i = i;
		this.j = j;
		this.k = k;
	}
	
	// Returns true if the job is contained in the set
	public boolean contains(int job)
	{
		if ( job >= i && job <= j )
		{
			if ( k==-1 ) // k=-1 means selecting all the elements in the range
				return true;
			else if ( Job.length(job) < Job.length(k) )
				return true;
			else if ( Job.length(job) == Job.length(k) && job<k )
				return true;
		}
		return false;
	}
	
	// Returns the number of jobs that are actually in the set
	public int count()
	{
		int n=0;
		for ( int i=this.i ; i<=this.j ; i++ )
			if ( this.contains(i) )
				n++;
		return n;
	}
	
	// Returns the job with the maximum processing time from the set
	// If multiple jobs fit the criterion, it returns the right-most in a
	// non-decreasing deadline ordering of the jobs.
	public int maxLength()
	{
		int rv = -22;
		for ( int i=this.i ; i<=this.j ; i++ )
			if ( this.contains(i) )
			{
				if ( rv==-22 )
					rv = i;
				else
					if ( Job.length(rv) <= Job.length(i))
						rv = i;
			}
		return rv;
	}
	
	// Returns the sum of the processing time of all the jobs contained in this set
	public int totalLength()
	{
		int total = 0;
		for ( int i=this.i ; i<=this.j ; i++ )
			if ( this.contains(i) )
				total += Job.length(i);
		return total;
	}
	
	// Returns the left-most
	// If the set is empty, a non-value is returned instead
	public int first()
	{
		for ( int i=this.i ; i<=this.j ; i++ )
			if ( contains(i) )
				return i;
		return -3;
	}
	
	public String toString()
	{
		return "S("+i+","+j+","+k+")";
	}
}


// Provides access to the data fro the jobs
class Job
{
	static int length(int job) { return algorithms.jobs[job][0]; }
	static int deadline(int job) { return algorithms.jobs[job][1]; }
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
		
		runAlgorithms();
	}
	
	// reads a problem, and outputs the result of both greedy and best-first
    public static void runAlgorithms ()
    {		
		// Run dynamic programming solution
//		System.out.print("Dynamic: ");
		System.out.println(Dynamic.search());
		
		// Run greedy solution
//		System.out.print("Greedy: ");
//		schedule s = greedy.greedy();
//		System.out.println(s.get_tardiness());
		
		// Run best first solution
//		boolean inTheMoodForCrashing = false;
//		if (inTheMoodForCrashing)
//		{
//			try {
//				System.out.print("Best first: ");
//				s = best_first_search.search();
//				System.out.println(s.get_tardiness());
//			} catch(Throwable e) {
//				// catches out of memory errors
//				//e.printStackTrace();
//				System.out.println("crashed");
//			}
//		}
	}
}
