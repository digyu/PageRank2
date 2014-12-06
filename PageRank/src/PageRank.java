//Page Rank - TCSS 554

import java.io.File;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Scanner;
import org.la4j.matrix.Matrix;
import org.la4j.matrix.sparse.CRSMatrix;
import org.la4j.vector.Vector;
import org.la4j.vector.dense.BasicVector;


/**
 * Is a page Rank Implementation that reads in a file that 
 * has an adjacency matrix with i j 1 format
 * 
 * @author Sanjeev Kamboj 
 * @date 12/4/14
 */
public class PageRank {
	
	/** Constant for matrix dimension*/
	private static final int NUM_ROW = 6;
	
	/** The damping factor*/
	private static final double DAMPING_FACTOR = .15;
	
	/** The name of the file to read */
	private static final String FILE_NAME = "Graph.txt";
	
	
	/** The main method*/
	public static void main(String[] args) throws FileNotFoundException {
		
		Scanner s = new Scanner(new File(FILE_NAME));
		
		Matrix matrix_A = makeMatrix(s);
		Matrix M = makeTransArray(matrix_A, s);
		runAlgorithm(M);
		
		
	}
	
	/**
	 * Reads in the file and makes the adjacency matrix
	 * 
	 * @param s - The scanner that reads in the file
	 * @return The adjacency matrix 
	 */
	public static Matrix makeMatrix(final Scanner s) {
		double[][] a = new double[NUM_ROW][NUM_ROW];
		
		while(s.hasNextLine()) {
			final String line_str = s.nextLine();
		    String[] split = line_str.split(" ");
		    
		    int row = Integer.parseInt(split[0]);
			int col = Integer.parseInt(split[1]);
			int bit = Integer.parseInt(split[2]);
			
			a[row-1][col-1] = bit;
				
		}
		final Matrix matrix_A = new CRSMatrix(a);
		
		return matrix_A;
		
	}
	
	/**
	 *  Makes the transition probability matrix from the
	 *  adjacency matrix
	 *  
	 * @param matrix_A
	 * @param s
	 * @return
	 */
	public static Matrix makeTransArray(final Matrix matrix_A, final Scanner s) {
		
		Matrix trans = new CRSMatrix(new double[NUM_ROW][NUM_ROW]);
		
		Queue<Integer> col = new LinkedList<Integer>();
		
		for(int i = 0; i < matrix_A.rows(); i++) {
			for(int j = 0; j < matrix_A.columns(); j++) {
				if(matrix_A.get(i, j) == 1.0) {
					col.add(j);
				}
			}
			trans = makeTransStepOne(i, col, trans);
			
		}
		
		//Step two, multiply by 1 - damping
		trans = trans.multiply((1.0-DAMPING_FACTOR));
		
		//Step three, add (damping/N)
		trans = trans.add((DAMPING_FACTOR / NUM_ROW));
		
		return trans;
	}
	
	/**
	 * If the row has no 1s, changes the value to 1/N
	 * Otherwise, changes the 1s to (1/number of ones in row)
	 * 
	 * @param i -The current row index
	 * @param col -A list of the column indices that contain 1
	 * @param trans The matrix to change
	 */
	private static Matrix makeTransStepOne(final int i, final Queue<Integer> col, 
			                             final Matrix trans) {
		final int num_ones = col.size();
		
		if(num_ones > 0) {
			while(!col.isEmpty()) {
				trans.set(col.remove(), i, (1.0 / num_ones));
			}
			
		} else {   //there are no ones in row
			for(int k = 0; k < trans.columns(); k++) {
				trans.set(k, i, (1.0 / (NUM_ROW)));
			} 
		}
		return trans;
	}
	
	/**
	 * Runs the Page Rank algorithm
	 * 
	 * @param m - the transition matrix with damping
	 */
	public static void runAlgorithm(final Matrix m) {
		Matrix p = m.copy();   // A in equation 
	    Vector prev = makeVectorV();  // V in equation - previous result
	    
		int t = 1;
		
	    Vector result = p.multiply(prev);  //A*v
	    
	    	
	    while(!equalVectors(result, prev)) {  //If old matrix and new one are same
	    	prev = result;
	    	result = p.multiply(result); 
	    	
	    	t++;
	    	
	    }
	    
	    double[] d_results = new double[result.length()];
	    
	    for(int i = 0; i < result.length(); i++) {
	    	d_results[i] = result.get(i);
	    	
	    }
	    
	    Arrays.sort(d_results);
	    printArray(d_results, t);
	   
		
	}
	
	/**
	 * Prints the Page Rank Array
	 * 
	 * @param array The array to print
	 * @param t The number of times the multiplication  (A^t * v)
	 */
	private static void printArray(final double[] array, final int t) {
		final DecimalFormat df = new DecimalFormat("#.###");
		
		System.out.println("At t = " + t +", Page Rank scores are: ");
		
		for(int i = array.length - 1; i >=0; i--) {
	        System.out.println(df.format(array[i]));
		}
	}
	
	/**
	 * Helper method to make the vector v where v = [1/N, 1N, ..]
	 * 
	 * @return THe Vector v
	 */
	private static Vector makeVectorV() {
		double[] v_array = new double[NUM_ROW];
		
		for(int i = 0; i < v_array.length ; i++) {
			v_array[i] = (1.0 /NUM_ROW);
			
		}
		Vector v = new BasicVector(v_array);
		
		return v;
	}
	
	/**
	 * Helper method to determine if two vectors are equal
	 * 
	 * @param v The vector to test against
	 * @param other The other vector to test again
	 * @return A boolean to see if the elements in the two vectors are the same
	 */
	private static boolean equalVectors(final Vector v, final Vector other) {
		boolean is_equal = true;
		
		if(v.length() == other.length()) {
			for(int i = 0; i < v.length(); i++) {
				if(v.get(i) != other.get(i)) {
					is_equal = false;
					break;
				}
			}
		} else {
			is_equal = false;
		}
		
		return is_equal;
		
	}
	
	
}
