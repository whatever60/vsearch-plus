package org.vsearchplus.rdp;

import edu.msu.cme.rdp.classifier.ClassificationResult;
import edu.msu.cme.rdp.classifier.GenusWordConditionalProb;
import edu.msu.cme.rdp.classifier.HierarchyTree;
import edu.msu.cme.rdp.classifier.RankAssignment;
import edu.msu.cme.rdp.classifier.ShortSequenceException;
import edu.msu.cme.rdp.classifier.TrainingInfo;
import edu.msu.cme.rdp.classifier.utils.ClassifierSequence;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.orientation.GoodWordIterator;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

/**
 * Native paired-end extension of the RDP Naive Bayes classifier.
 *
 * This class keeps the stock RDP scoring mechanics and extends the query unit
 * from one sequence to one ordered pair of anchor sequences.
 */
public final class PairedNaiveBayesClassifier {
  private static final int NUM_OF_RUNS = 100;
  private static final int MAX_NUM_OF_WORDS = 5000;

  private final TrainingInfo trainingInfo;
  private final long bootstrapSeed;
  private final Random randomGenerator;
  private final Random randomSelectGenera;

  private float[][] querySeqWordProbArr;
  private float[] accumulateProbArr;

  /**
   * Create a paired classifier from an already-loaded stock training info
   * object.
   */
  public PairedNaiveBayesClassifier(final TrainingInfo trainingInfo,
                                    final long bootstrapSeed) {
    this.trainingInfo = trainingInfo;
    this.bootstrapSeed = bootstrapSeed;
    final int nodeListSize = trainingInfo.getGenusNodeListSize();
    this.querySeqWordProbArr = new float[MAX_NUM_OF_WORDS][nodeListSize];
    this.accumulateProbArr = new float[nodeListSize];
    this.randomGenerator = new Random(bootstrapSeed);
    this.randomSelectGenera = new Random(bootstrapSeed);
  }

  /**
   * Classify one paired query as a single TAV object.
   */
  public ClassificationResult
  classifyPair(final String tavSequenceId, final ClassifierSequence leftInput,
               final ClassifierSequence rightInput, final int minBootstrapWords)
      throws IOException {
    final boolean leftWasReversed = trainingInfo.isSeqReversed(
        leftInput.getWordIndexArr(), leftInput.getGoodWordCount());
    final ClassifierSequence left =
        leftWasReversed ? leftInput.getReversedSeq() : leftInput;

    final boolean rightWasReversed = trainingInfo.isSeqReversed(
        rightInput.getWordIndexArr(), rightInput.getGoodWordCount());
    final ClassifierSequence right =
        rightWasReversed ? rightInput.getReversedSeq() : rightInput;

    final int[] leftWordIndexArr = left.getWordIndexArr();
    final int[] rightWordIndexArr = right.getWordIndexArr();
    final int leftGoodWordCount = left.getGoodWordCount();
    final int rightGoodWordCount = right.getGoodWordCount();
    final int totalGoodWordCount = leftGoodWordCount + rightGoodWordCount;

    final int totalLength =
        left.getSeqString().length() + right.getSeqString().length();
    if (totalLength < edu.msu.cme.rdp.classifier.Classifier.MIN_SEQ_LEN) {
      throw new ShortSequenceException(
          tavSequenceId,
          "ShortSequenceException: Paired sequence with recordID=" +
              tavSequenceId + " is less than " +
              edu.msu.cme.rdp.classifier.Classifier.MIN_SEQ_LEN);
    }

    if (totalGoodWordCount <
        edu.msu.cme.rdp.classifier.Classifier.MIN_GOOD_WORDS) {
      throw new ShortSequenceException(
          tavSequenceId,
          "ShortSequenceException: Paired sequence with recordID=" +
              tavSequenceId + " does not have enough valid words, minimum " +
              edu.msu.cme.rdp.classifier.Classifier.MIN_GOOD_WORDS +
              " are required");
    }

    ensureCapacity(totalGoodWordCount, trainingInfo.getGenusNodeListSize());

    populateWordProbabilityRows(leftWordIndexArr, leftGoodWordCount, 0);
    populateWordProbabilityRows(rightWordIndexArr, rightGoodWordCount,
                                leftGoodWordCount);

    final int determinedNodeIndex =
        determineBestNodeFromAllWords(totalGoodWordCount);
    final HashMap<HierarchyTree, Float> confidenceCountMap =
        initConfidenceMap(determinedNodeIndex);

    final int numSelections = Math.max(
        totalGoodWordCount / GoodWordIterator.getWordsize(), minBootstrapWords);

    randomGenerator.setSeed(bootstrapSeed);
    randomSelectGenera.setSeed(bootstrapSeed);
    runBootstrap(totalGoodWordCount, numSelections, confidenceCountMap);

    final List<RankAssignment> finalAssignments =
        finalizeAssignments(determinedNodeIndex, confidenceCountMap);

    final Sequence querySequence =
        new Sequence(tavSequenceId, "", left.getSeqString());
    return new ClassificationResult(querySequence, leftWasReversed,
                                    finalAssignments,
                                    trainingInfo.getHierarchyInfo());
  }

  /**
   * Grow temporary probability buffers when paired word count exceeds defaults.
   */
  private void ensureCapacity(final int wordCount, final int nodeListSize) {
    if (wordCount > querySeqWordProbArr.length) {
      querySeqWordProbArr = new float[wordCount][nodeListSize];
    }
  }

  /**
   * Populate classifier word-probability rows for one anchor.
   */
  private void populateWordProbabilityRows(final int[] wordIndexArr,
                                           final int goodWordCount,
                                           final int rowOffset) {
    final int nodeListSize = trainingInfo.getGenusNodeListSize();

    for (int offset = 0; offset < goodWordCount; offset++) {
      final int wordIndex = wordIndexArr[offset];
      final float wordPrior = trainingInfo.getLogWordPrior(wordIndex);
      final int outRow = rowOffset + offset;

      for (int node = 0; node < nodeListSize; node++) {
        querySeqWordProbArr[outRow][node] =
            wordPrior - trainingInfo.getLogLeaveCount(node);
      }

      final int start = trainingInfo.getStartIndex(wordIndex);
      final int stop = trainingInfo.getStopIndex(wordIndex);
      for (int n = start; n < stop; n++) {
        final GenusWordConditionalProb prob =
            trainingInfo.getWordConditionalProbObject(n);
        querySeqWordProbArr[outRow][prob.getGenusIndex()] =
            prob.getProbability();
      }
    }
  }

  /**
   * Determine the top genus node from full paired evidence.
   */
  private int determineBestNodeFromAllWords(final int totalGoodWordCount) {
    final int nodeListSize = trainingInfo.getGenusNodeListSize();

    for (int row = 0; row < totalGoodWordCount; row++) {
      for (int node = 0; node < nodeListSize; node++) {
        accumulateProbArr[node] += querySeqWordProbArr[row][node];
      }
    }

    float maxPosteriorProb = Float.NEGATIVE_INFINITY;
    int determinedNodeIndex = 0;
    for (int node = 0; node < nodeListSize; node++) {
      if (accumulateProbArr[node] > maxPosteriorProb) {
        determinedNodeIndex = node;
        maxPosteriorProb = accumulateProbArr[node];
      }
      accumulateProbArr[node] = 0.0f;
    }
    return determinedNodeIndex;
  }

  /**
   * Initialize confidence counters on the determined lineage path.
   */
  private HashMap<HierarchyTree, Float>
  initConfidenceMap(final int determinedNodeIndex) {
    final HashMap<HierarchyTree, Float> confidenceCountMap =
        new HashMap<HierarchyTree, Float>();
    HierarchyTree node = trainingInfo.getGenusNodebyIndex(determinedNodeIndex);
    while (node != null) {
      confidenceCountMap.put(node, 0.0f);
      node = node.getParent();
    }
    return confidenceCountMap;
  }

  /**
   * Run paired bootstrap on the combined R1+R2 word pool.
   */
  private void
  runBootstrap(final int totalGoodWordCount, final int numSelections,
               final HashMap<HierarchyTree, Float> confidenceCountMap) {
    final int nodeListSize = trainingInfo.getGenusNodeListSize();

    for (int run = 0; run < NUM_OF_RUNS; run++) {
      for (int draw = 0; draw < numSelections; draw++) {
        final int randomIndex = randomGenerator.nextInt(totalGoodWordCount);
        for (int node = 0; node < nodeListSize; node++) {
          accumulateProbArr[node] += querySeqWordProbArr[randomIndex][node];
        }
      }

      float maxPosteriorProb = Float.NEGATIVE_INFINITY;
      int bestNodeIndex = 0;
      boolean tied = false;

      for (int node = 0; node < nodeListSize; node++) {
        if (accumulateProbArr[node] >= maxPosteriorProb) {
          if (accumulateProbArr[node] > maxPosteriorProb) {
            bestNodeIndex = node;
            maxPosteriorProb = accumulateProbArr[node];
            tied = false;
          } else {
            tied = true;
          }
        }
      }

      if (tied) {
        final ArrayList<Integer> possibleSet = new ArrayList<Integer>();
        for (int node = 0; node < nodeListSize; node++) {
          if (accumulateProbArr[node] == maxPosteriorProb) {
            possibleSet.add(node);
          }
        }
        bestNodeIndex =
            possibleSet.get(randomSelectGenera.nextInt(possibleSet.size()));
      }

      for (int node = 0; node < nodeListSize; node++) {
        accumulateProbArr[node] = 0.0f;
      }

      incrementPathConfidence(trainingInfo.getGenusNodebyIndex(bestNodeIndex),
                              confidenceCountMap);
    }
  }

  /**
   * Add one bootstrap vote to the lineage path if it is part of determined
   * lineage.
   */
  private void incrementPathConfidence(
      final HierarchyTree nodeStart,
      final HashMap<HierarchyTree, Float> confidenceCountMap) {
    HierarchyTree node = nodeStart;
    while (node != null) {
      final Float previous = confidenceCountMap.get(node);
      if (previous != null) {
        confidenceCountMap.put(node, previous + 1.0f);
      }
      node = node.getParent();
    }
  }

  /**
   * Convert bootstrap counts to rank assignments in root-to-leaf order.
   */
  private List<RankAssignment>
  finalizeAssignments(final int determinedNodeIndex,
                      final HashMap<HierarchyTree, Float> confidenceCountMap) {
    final List<RankAssignment> finalAssignments =
        new ArrayList<RankAssignment>();
    HierarchyTree node = trainingInfo.getGenusNodebyIndex(determinedNodeIndex);
    while (node != null) {
      final float confidence =
          confidenceCountMap.get(node) / (float)NUM_OF_RUNS;
      finalAssignments.add(0, new RankAssignment(node, confidence));
      node = node.getParent();
    }
    return finalAssignments;
  }
}
