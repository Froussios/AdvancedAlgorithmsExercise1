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
	
	public static int search()
	{
		Problem problem = Problem.LoadProblem();
		
		int n = problem.size();
		memory = new int[n][n][n][n];
		for (int i1=0 ; i1<n ; i1++ )
			for (int i2 = 0; i2 < n; i2++)
				for (int i3 = 0; i3 < n; i3++)
					for (int i4 = 0; i4 < n; i4++)
						memory[i1][i2][i3][i4] = -1;
		
		ArrayList<Integer> sequence = new ArrayList<Integer>();
		int tardiness = recursive(new Subset(0, algorithms.num_jobs-1, -1), 0, 0, sequence);
		
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
		
		return tardiness;
	}
	
	static int run(Problem problem, int t)
	{
		return -1;
	}
	
	static int hits = 0;
	static int misses = 0;
	static int recursive(Subset subset, int t, int depth, ArrayList<Integer> rvSeq)
	{
		// Count items
		int n=0;
		for ( int i=subset.i ; i<=subset.j ; i++ )
			if ( subset.contains(i) )
				n++;
		
		// Subset is empty
		if (n == 0)
			return 0;
		
		// Find k : k belongs to subset and has maximum length
		int k = subset.maxLength();
		
		int best = Integer.MAX_VALUE;
		ArrayList<Integer> best_sequence = null;
		
		for ( int split=k ; split <= subset.j ; split++ )
		{
			Subset subset1 = new Subset(subset.i,split,k);
			Subset subset2 = new Subset(split+1,subset.j,k);
			
			ArrayList<Integer> seqPre = new ArrayList<Integer>();
			ArrayList<Integer> seqPost = new ArrayList<Integer>();
			
			int completion_k = t + subset1.totalLength() + Job.length(k);
			
			int tard1 = recursive(subset1, t, depth++, seqPre);
			int tard2 = recursive(subset2, completion_k, depth++, seqPost);
			int tardk = Math.max(0, completion_k - Job.deadline(k));
			
			int tard = tard1 + tardk + tard2;
			
			if ( tard < best )
			{
				best = tard;
				best_sequence = seqPre;
				best_sequence.add(k);
				best_sequence.addAll(seqPost);
			}
		}
		
		if ( best == Integer.MAX_VALUE )
			System.out.println("Failed to find solution: " + subset + " k:" + k);
		
		rvSeq.addAll(best_sequence);
		
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
		ArrayList<Integer> seqBest = null;
		
		for ( int d=0 ; d <= n-k-1 ; d++ )
		{
			Problem subproblem1 = new Problem(k+d/*+1*/);
			for ( int i=0 ; i<k+d+1 ; i++ )
				if ( i!=k )
					subproblem1.add(problem.get(i));
			
			Problem subproblem2 = new Problem(n - (k+d));
			for ( int i=k+d+1 ; i<n ; i++ )
				subproblem2.add(problem.get(i));
			
			ArrayList<Integer> seqPre = new ArrayList<Integer>();
			ArrayList<Integer> seqPost = new ArrayList<Integer>();
			
			int tard1 = recursive(subproblem1, start, seqPre);
			int finish_k = start + subproblem1.totalLength() + problem.get(k).length;
			int tard2 = recursive(subproblem2, finish_k, seqPost);
			
			int tard_k = tard1 + tard2 + Math.max(finish_k-problem.get(k).deadline, 0);
			
			if (tard_k < best)
			{
				best = tard_k;
				seqBest = seqPre;
				seqBest.add(problem.get(k).job);
				seqBest.addAll(seqPost);
			}
		}
		
		//System.out.println("\tJob " + problem.get(k) + " start: " + best_preceeding.totalLength() + " tardiness " + Math.max(0, start - problem.get(k).deadline));
		
		rvSeq.addAll(seqBest);
		
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

class Subset
{
	public int i,j,k;
	public Subset(int start/*inclusive*/, int end/*inclusive*/, int lessthan)
	{
		assert(start >= 0);
		assert(end >= start);
		this.i = start;
		this.j = end;
		this.k = lessthan;
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
		return "" + job + " (" +length+","+deadline+")";
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

	// reads a problem, and outputs the result of both greedy and best-first
    public static void main (String args[]) {
		read_problem(args[0]);
		
		// Sort problem in non-decreasing order on deadline
		Arrays.sort(algorithms.jobs, new Comparator<int[]>() {
				public int compare(int[] f1, int[] f2)
					{
		            	return f1[1] - f2[1];
					}
			});
		
		
		{
			// DEBUGGING: Print problem
			System.out.println("Problem: (length,due)");
			for (int[] i : algorithms.jobs )
				System.out.print( "(" + i[0] + "," +i [1] + "), ");
			System.out.println();
		}
		
		// Run brute force version of dynamic solution
		System.out.print("Brute force: ");
		System.out.println(BruteForce.search());
		
		// Run recursive version of dynamic solution
		System.out.print("Recursive: ");
		System.out.println(Recursive.search());
		
		// Run dynamic programming solution
		System.out.print("Dynamic: ");
		System.out.println(Dynamic.search());
		
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
