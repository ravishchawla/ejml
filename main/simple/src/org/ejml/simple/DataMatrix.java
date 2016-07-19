package org.ejml.simple;

import org.ejml.simple.SimpleMatrix;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * A Multi-Purpose Data-Oriented Matrix class that is useful for storing data corresponding to 
 * name-based columns. DataMatrix implements Iterable, and can be used to iterate
 * over in a sliding or a fixed-window of arbitrary size. It also implements
 * many important linear algebra functions that are useful for computation in matrices.
 * Created by Ravish Chawla on 6/28/2016.
 */
public class DataMatrix extends SimpleMatrix implements Iterable<DataMatrix> {


	/**
	 * Default value of a Double used when provided values are invalid
	 */
	public static final double DEFAULT_DOUBLE = 0.0;

	/**
	* Default value of an Integer used when provided values are invalid
	*/
	public static final int DEFAULT_INT = -1;

	/**
	* Default value of a String used whe provided values are invalid
	*/
	public static final String DEFAULT_STRING = "DEFAULT_STRING";

	/**
	* Used for stroring names of columns in the data
	*/
	private String[] columns;

	/**
	* Iterator class implemetns Iterator. Used for iterating over the data in a 
	* window of data in slding of fixed chunks of data.
	*/

	private class Iterator implements java.util.Iterator<DataMatrix> {

		/**
		 * Used for holding the position where the iterator is positioned at.
		 */
		private int cursor;

		/**
		 * Holds the length of the datamatrix.
		 */
		private int length;

		/**
		 * Holds the size of the window over which to iterate.
		 */
		private int windowSize;

		/**
		 * The number of columns in the datamatrix.
		 */
		private int numCols;

		/**
		 * The amount by which to increment at each next call.
		 */
		private int increment;

		/**
		 * Constructor to instantiate Iterator with windowsize and the option to use sliding or fixed iteration.
		 * @param windowSize The size of the window over which to iterate. The smaller of the WindowSize and Size of the dataset is used.
		 * @param sliding If true, use a sliding window. If false, use fixed.
		 */
		public Iterator(int windowSize, boolean sliding) {
			this.cursor = 0;

			DataMatrix mm = DataMatrix.this;
			this.length  = mm.numRows();
			this.numCols = mm.numCols();
			this.windowSize = Math.min(windowSize, length);
			this.increment = (sliding ? 1 : windowSize);
		}

		/**
		 * Move the iterator to the next section of the data of provided WindowSize.
		 * @return Next iterable section of the DataMatrix. Null if the iterator has reached the end.
		 */
		public DataMatrix next() {

			if(!this.hasNext()) {
				return null;
			}

			DataMatrix ret = DataMatrix.this.extractMatrix(this.cursor, this.cursor + this.windowSize, 0, this.numCols);
			this.cursor = this.cursor + this.increment;
			return ret;
		}

		/**
		 * Determine if the iterator has reached the end of the DataMatrix.
		 * @return True if there are more iterable sections. False otherwise.
		 */
		public boolean hasNext() {
			return (this.cursor <= (this.length - this.windowSize));
		}

		/**
		 * This implementation of the iterator does not implement remove()
		 */
		public void remove() {

		}
	}

	/**
	 * Create a DataMatrix of fixed size.
	 * @param numRows The number of rows.
	 * @param numCols The number of columns.
	 */
	public DataMatrix(int numRows, int numCols) {
		super(numRows, numCols);
		this.columns = new String[numCols];
	}

	/**
	 * Create a DataMatrix Column Vector with provided values.
	 * @param values The values from which to instantiate the DataMatrix.
	 */
	public DataMatrix(Float[] values) {
		super(values.length, 1);
		for(int i = 0; i < values.length; i++) {
			this.set(i, 0, values[i]);
		}

		this.columns = new String[0];
	}

	/**
	 * Create a DataMatrix from an existing SimpleMatrix.
	 * @param m The SimpleMatrix from which to instantiate the DataMatrix.
	 */
	protected DataMatrix(SimpleMatrix m) {
		this.mat = m.getMatrix();
		this.columns = null;
	}

	/**
	 * Create a DataMatrix from an existing DenseMatrix64F.
	 * @param m The DenseMatrix64F from which to instantiate the DataMatrix.
	 */
	protected DataMatrix(DenseMatrix64F m) {
		this.mat = m;
		this.columns = null;
	}

	/**
	 * Create an Iterator of default values.
	 * @return An Iterator of WindowSize 10 and Sliding Windows.
	 */
	public java.util.Iterator<DataMatrix> iterator() {
		return new Iterator(10, true);
	}

	/**
	 * Create an Iterator with given Window Size and Sliding Window option.
	 * @param windowSize The Window Size
	 * @param sliding Option. If True, use Sliding Windows. Fixed otherwise.
	 * @return An Iterator with the provided options.
	 */
	public java.util.Iterator<DataMatrix> iterator(int windowSize, boolean sliding) {
		return new Iterator(windowSize, sliding);
	}

	/**
	 * Create an Iterator with given Window Size and default Sliding Window option.
	 * @param windowSize The Window Size
	 * @return An Iterator with the provided WindowSize and a Sliding Window.
	 */
	public java.util.Iterator<DataMatrix> iterator(int windowSize) {
		return new Iterator(windowSize, true);
	}

	/**
	 * Create an Iterator with given Sliding Window option and default Window Size
	 * @param sliding If true, use Sliding Windows. Fixed otherwisel
	 * @return An Iterator with provided Sliding Window option and a Window Size of 10.
	 */
	public java.util.Iterator<DataMatrix> iterator(boolean sliding) {
		return new Iterator(10, sliding);
	}

	/**
	 * Calculate the mean of the entire DataMatrix over all cells.
	 * @return The Mean of the DataMatrix.
	 */
	public double mean() {
		double sum = 0;
		for(int i = 0; i < this.numRows(); i++) {
			for(int j = 0; j < this.numCols(); j++) {
				sum = sum + this.get(i, j);
			}
		}

		return (sum / this.getNumElements());
	}

	/**
	 * Calculate the mean of all cells in a single column in the DataMatrix.
	 * @param column The Column for which to calculate the mean for.
	 * @return The Mean of the column in the DataMatrix. DEFAULT_DOUBLE if column is not in the matrix.
	 */
	public double mean(String column) {
		double sum = 0;
		int columnIndex = this.getColumnIndex(column);

		if(columnIndex == DataMatrix.DEFAULT_INT) {
			return DataMatrix.DEFAULT_DOUBLE;
		}

		for(int i = 0; i < this.numRows(); i++) {
			sum = sum + this.get(i, columnIndex);
		}

		return (sum / this.getNumElements());
	}

	/**
	 * Calculate the Standard Deviation of all cells in the DataMatrix.
	 * @return The Standard Deviation of the DataMatrix.
	 */
	public double std() {
		double avg = this.mean();
		double squaredSum = 0;

		int size = this.getNumElements();
		if(size <= 1) {
			return 0;
		}

		for(int i = 0; i < size; i++) {
			double val = this.get(i);
			squaredSum = ((val - avg) * (val - avg));
		}

		squaredSum = squaredSum / (size - 1);
		return Math.sqrt(squaredSum);
	}

	/**
	 * Calculate a differentating Matrix, where the value of each row in the Differentating Matrix is
	 * the difference of value in that row minus the value in the nth previous row.
	 * @param row Number of rows by which to calculate the difference by.
	 * @return A Differentiating Matrix of size (NumRows - row)
	 */
	public DataMatrix diffByRowMatrix(int row) {
		DataMatrix dataMatrix = new DataMatrix(this.extractMatrix(row, this.numRows(), 0, this.numCols()).minus(
				this.extractMatrix(0, this.numRows() - row, 0, this.numCols())));

		dataMatrix.columns = this.columns;

		return dataMatrix;
	}

	/**
	 * Extract a part of the DataMatrix. The columns of the new matrix are a subset of the original DataMatrix within the new boundaries.
	 * @param y0 The starting Row
	 * @param y1 The ending Row
	 * @param x0 The starting Column
	 * @param x1 The ending Column
	 * @return A DataMatrix that contains elements within the provided points.
	 */
	@Override
	public DataMatrix extractMatrix(int y0, int y1, int x0, int x1) {
		if(this.numRows() == 0 || this.numCols() == 0) {
			System.out.println("DataMatrix is of form 0xn or nx0");
			return this;
		}
		DataMatrix dataMatrix = new DataMatrix(super.extractMatrix(y0, y1, x0, x1));
		dataMatrix.columns = new String[dataMatrix.numCols()];

		if(this.columns == null) {
			System.out.println("DM Column is NULL");
			return dataMatrix;
		}

		for(int i = x0; i < x1; i++) {
			dataMatrix.columns[i - x0] = this.columns[i];
		}

		return dataMatrix;
	}

	/**
	 * Extract Only specific columns in the DataSet as a new DataMatrix.
	 * @param columns The columns to extract.
	 * @return The DataMatrix of the corresponding columns.
	 */
	public DataMatrix extractColumns(String... columns) {
		if(columns.length == 0) {
			return null;
		}

		int[] columnIndices = new int[columns.length];
		for(int i = 0; i < columns.length; i++) {
			columnIndices[i] = this.getColumnIndex(columns[i]);
			if(columnIndices[i] == DataMatrix.DEFAULT_INT) {
				System.out.println("Provided Column is not in the DataMatrix");
				return this;
			}
		}

		DataMatrix first = this.extractMatrix(0, this.numRows(), columnIndices[0], columnIndices[0] + 1);

		for(int i = 1; i < columnIndices.length; i++) {
			first = first.union(this.extractMatrix(0, this.numRows(), columnIndices[i], columnIndices[i] + 1));
		}

		first.setColumns(columns);

		return first;
	}

	/**
	 * Calculates the Transpose of the DataMatrix.
	 * The Columns of the Transposed matrix are going to be the same as the original columns.
	 * @return Transposed matrix.
	 */
	public DataMatrix transpose() {
		DataMatrix transpose = new DataMatrix(super.transpose());
		transpose.columns = this.columns;
		return transpose;
	}

	/**
	 * Fill all cells in the DataMatrix with a scalar value.
	 * @param value The value that should be filled.
	 */
	public void fill(double value) {
		for(int i = 0; i < this.numRows(); i++) {
			for(int j = 0; j < this.numCols(); j++) {
				this.set(i, j, value);
			}
		}
	}

	/**
	 * Scalar Power operation on all cells.
	 * @param pow The Power value to use.
	 * @return A DataMatrix with all elements raised to the power.
	 */
	@Override
	public DataMatrix elementPower(double pow) {
		DataMatrix dataMatrix = new DataMatrix(super.elementPower(pow));
		dataMatrix.columns = this.columns;
		return dataMatrix;
	}

	/**
	 * Add Two Matrices.
	 * @param to The DataMatrix to add.
	 * @return The DataMatrix sum.
	 */
	@Override
	public DataMatrix plus(SimpleMatrix to) {
		DataMatrix dataMatrix = new DataMatrix(super.plus(to));
		dataMatrix.columns = this.columns;
		return dataMatrix;
	}

	/**
	 * Multiply Two Matrices using element-wise multiplication.
	 * @param with The Matrix to multiply.
	 * @return The DataMatrix product.
	 */
	@Override
	public DataMatrix mult(SimpleMatrix with) {
		DataMatrix dataMatrix = new DataMatrix(super.mult(with));
		dataMatrix.columns = this.columns;
		return dataMatrix;
	}

	/**
	 * Cross Product of DataMatrix Triplet Vectors.
	 * @param with The DataMatrix by which to cross.
	 * @return The Cross Product.
	 */
	public DataMatrix tripletVectorCrossProduct(DataMatrix with) {
		if(!this.isTriplet()) {
			System.out.println("Library only supports Cross Product of three-element vectors");
			return null;
		}

		DataMatrix resultant = new DataMatrix(1, 3);
		DataMatrix first = this.getStandardTripletForm();
		DataMatrix second = with.getStandardTripletForm();

		resultant.set(0, 0, first.get(0, 1) * second.get(0, 2) - first.get(0, 2) * second.get(0, 1));
		resultant.set(0, 1, first.get(0, 2) * second.get(0, 0) - first.get(0, 0) * second.get(0, 2));
		resultant.set(0, 2, first.get(0, 0) * second.get(0, 1) - first.get(0, 1) * second.get(0, 0));

		return resultant;

	}

	/**
	 * Calculate the Sigmoid Function on the elements of the DataMatrix.
	 * @return DataMatrix where each element has sigmoid applied.
	 */
	public DataMatrix sigmoid() {

		//Ha!
		SimpleMatrix s = this.negative().elementExp().plus(1).elementPower(-1);
		DataMatrix dataMatrix = new DataMatrix(s);
		dataMatrix.columns = this.columns;
		return dataMatrix;
	}

	/**
	 * Calculate a Binary Split on a provided value.
	 * @param splitValue The value by which to split by.
	 * @return A DataMatrix such that all values that are less than splitValue are 0, all values greater are 1.
	 */
	public DataMatrix binarySplit(double splitValue) {
		SimpleMatrix dataMatrix = this.copy();

		for(int i = 0; i < dataMatrix.numRows(); i++) {
			for(int j = 0; j < dataMatrix.numCols(); j++) {
				dataMatrix.set(i, j, dataMatrix.get(i, j) < splitValue ? 0 : 1);
			}
		}

		DataMatrix returnMatrix = new DataMatrix(dataMatrix);
		returnMatrix.columns = this.columns;

		return returnMatrix;
	}

	/**
	 * Majority element of the Matrix.
	 * @return the element that apears most often in the Matrix.
	 */
	public double getMajorityElement() {
		Map<Double, Integer> counts = new HashMap<Double, Integer>();
		Double majority = this.get(0, 0);

		for(int i = 0; i < this.numRows(); i++) {
			for(int j = 0; j < this.numCols(); j++) {
				double key = this.get(i, j);

				if(!counts.containsKey(key)) {
					counts.put(key, 1);
				} else {
					counts.put(key, counts.get(key) + 1);
				}

				if(counts.get(key) > majority) {
					majority = key;
				}
			}
		}

		return  majority;
	}

	/**
	 * The Minimum element in the DataMatrix.
	 * @return Smallest element.
	 */
	public double getMinElement() {
		double min = Integer.MAX_VALUE;
		for(int row = 0; row < this.numRows(); row++) {
			for(int col = 0; col < this.numCols(); col++) {
				if(this.get(row, col) < min) {
					min = this.get(row, col);
				}
			}
		}

		if(min == Integer.MAX_VALUE) {
			return DataMatrix.DEFAULT_DOUBLE;
		} else {
			return min;
		}
	}

	/**
	 * The Maximum element in the DataMatrix.
	 * @return Larget Element.
	 */
	public double getMaxElement() {
		double max = Integer.MIN_VALUE;
		for(int row = 0; row < this.numRows(); row++) {
			for(int col = 0; col < this.numCols(); col++) {
				if(this.get(row, col) > max) {
					max = this.get(row, col);
				}
			}
		}

		if(max == Integer.MIN_VALUE) {
			return DataMatrix.DEFAULT_DOUBLE;
		} else {
			return max;
		}
	}

	/**
	 * The Range of the DataMatrix.
	 * @return The Difference between the largest and smalleset element in the Matrix.
	 */
	public double getRange() {
		double min = Integer.MIN_VALUE;
		double max = Integer.MAX_VALUE;
		for(int row = 0; row < this.numRows(); row++) {
			for(int col = 0; col < this.numCols(); col++) {
				if(this.get(row, col) > max) {
					max = this.get(row, col);
				}

				if(this.get(row, col) < min) {
					min = this.get(row, col);
				}
			}
		}

		if(max == Integer.MIN_VALUE || min == Integer.MAX_VALUE) {
			return DataMatrix.DEFAULT_DOUBLE;
		} else {
			return (max - min);
		}
	}

	/**
	 * The sum of all elements in the datamatrix.
	 * @return The sum.
	 */
	public double sum() {
		double sum = 0;
		for(int row = 0; row < this.numRows(); row++) {
			for(int col = 0; col < this.numCols(); col++) {
				sum = sum + this.get(row, col);
			}
		}
		return sum;
	}

	/**
	 * Add a unit colummn to the DataMatrix.
	 * @return A DataMatrix of (nx1 by m) size, where the first column values are all set to 1.
	 */
	public DataMatrix addUnitColumn() {
		DataMatrix ones = new DataMatrix(this.numRows(), 1);
		ones.fill(1);

		return ones.union(this);
	}

	/**
	 * The Union of Two DataMatrices. Both Matrices must have the same number of rows.
	 * @param second The DataMatrix to union with.
	 * @return A Unified DataMatrix. The columns of the new DataMatrix are consecutively combined as well.
	 */
	public DataMatrix union(DataMatrix second) {

		if(second == null) {
			return this;
		}

		if(this.numRows() != second.numRows()) {
			System.out.println("Number of rows mismatch in Union method");
			return null;
		}

		boolean setColumns = true;
		DataMatrix dataMatrix = new DataMatrix(this.numRows(), this.numCols() + second.numCols());
		for(int row = 0; row < this.numRows(); row++) {
			for(int col = 0; col < this.numCols(); col++) {
				dataMatrix.set(row, col, this.get(row, col));
				if(setColumns) {
					dataMatrix.columns[col] = this.columns[col];
				}
			}

			for(int col = this.numCols(); col < this.numCols() + second.numCols(); col++) {
				dataMatrix.set(row, col, second.get(row, col - this.numCols()));
				if(setColumns) {
					dataMatrix.columns[col] = second.columns[col - this.numCols()];
				}
			}

			setColumns = false;
		}

		return dataMatrix;
	}

	/**
	 * The Union of more than one DataMatrix. All DataMatrices must have the same number of
	 * rows.
	 * @param matrices The Matrices to unify with.
	 * @return The Unified DataMatrix. The columns of the new DataMatrix are also
	 * consecutively combined.
	 */
	public DataMatrix union(List<DataMatrix> matrices) {
		if(matrices.size() < 0) {
			return this;
		}

		DataMatrix finalUnion = this;
		for(DataMatrix matrix : matrices) {
			finalUnion = finalUnion.union(matrix);
			if(finalUnion == null) {
				return null;
			}
		}

		return finalUnion;
	}

	/**
	 * Concatination of two DataMatrices. The number of columns of both matrices must be the same.
	 * @param second The DataMatrix to concat with.
	 * @return The Concatinated DataMatrix.
	 */
	public DataMatrix concat(DataMatrix second) {

		if(second == null) {
			return this;
		}

		if(this.numCols() != second.numCols()) {
			System.out.println("Number of cols mismatch in Concat method");
			return null;
		}

		DataMatrix dataMatrix = new DataMatrix(this.numRows() + second.numRows(), this.numCols());

		boolean setColumns = true;
		for(int row = 0; row < this.numRows(); row++) {
			for(int col = 0; col < this.numCols(); col++) {
				dataMatrix.set(row, col, this.get(row, col));
				if(setColumns) {
					dataMatrix.columns[col] = this.columns[col];
				}
			}
			setColumns = false;
		}

		for(int row = this.numRows(); row < this.numRows() + second.numRows(); row++) {
			for(int col = 0; col < this.numCols(); col++) {
				dataMatrix.set(row, col, second.get(row - this.numRows(), col));
			}
		}

		return dataMatrix;
	}

	/**
	 * Concatination of multiple DataMatrices. The number of columns of all matrices must be the same.
	 * @param matrices The DataMatrix to concat.
	 * @return The concatinated DataMatrix.
	 */
	public DataMatrix concat(List<DataMatrix> matrices) {
		if(matrices.size() < 0) {
			return this;
		}

		DataMatrix finalConcat = this;
		for(DataMatrix matrix : matrices) {
			finalConcat = finalConcat.concat(matrix);
			if(finalConcat == null) {
				return null;
			}
		}

		return finalConcat;
	}

	/**
	 * Get the Index of a specified column name.
	 * @param columnName The column for which to get the index of.
	 * @return The column index.
	 */
	public int getColumnIndex(String columnName) {
		for(int i = 0; i < this.columns.length; i++) {
			if(this.columns[i].equals(columnName)) {
				return i;
			}
		}
		return DataMatrix.DEFAULT_INT;
	}

	/**
	 * Get the Column name at a specified index.
 	 * @param column The index for which to get the column name of.
	 * @return The column name.
	 */
	public String getColumn(int column) {
		if(column > this.numCols()) {
			return DataMatrix.DEFAULT_STRING;
		}

		return this.columns[column];
	}

	/**
	 * Get the column names as a String Array.
	 * @return
	 */
	public String[] getColumns() {
		return this.columns;
	}

	/**
	 * Set the columns of the DataMatrix. The number of columns must be the same as the size of the array.
	 * @param columns
	 */
	public void setColumns(String[] columns) {
		if (columns.length != this.numCols()) {
			System.out.println("Number of index columns provided do not match matrix columns");
			return;
		}

		this.columns = columns;
	}

	/**
	 * Get the Shape of the DataMatrix in the format (#rows, #columns)
	 * @return
	 */
	public int[] shape() {
		int[] point = {this.numRows(), this.numCols()};
		return point;
	}

	/**
	 * Check whether the DataMatrix is of Vector format.
	 * @return True if vector. False otherwise.
	 */
	public boolean isVector() {
		int[] point = this.shape();
		return (point[0] ==1 || point[1] == 1);
	}

	/**
	 * Check whether the DataMatrix is of Triplet Vector Format.
	 * @return True if Triplet. False otherwise.
	 */
	public boolean isTriplet() {
		int[] point = this.shape();
		return (this.isVector() && (point[0] == 3 || point[1] == 3));
	}

	/**
	 * Get a Triplet Vector in standard Row-Vector format.
	 * @return DataMatrix Vector in standard format.
	 */
	public DataMatrix getStandardTripletForm() {
		if(!isTriplet()) {
			System.out.println("Matrix not in triplet format (3 x 1 or 1 x 3)");
			return null;
		}

		if(this.numRows() == 1) {
			return this;
		} else {
			return this.transpose();
		}
	}
}