package org.vsearchplus.rdp;

import edu.msu.cme.rdp.classifier.ClassificationResult;
import edu.msu.cme.rdp.classifier.ShortSequenceException;
import edu.msu.cme.rdp.classifier.TrainingInfo;
import edu.msu.cme.rdp.classifier.io.ClassificationResultFormatter;
import edu.msu.cme.rdp.classifier.utils.ClassifierFactory;
import edu.msu.cme.rdp.classifier.utils.ClassifierSequence;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.HashMap;

/**
 * CLI entrypoint for native paired-end TAV taxonomy assignment using stock RDP internals.
 */
public final class PairedClassifierMain {
    private PairedClassifierMain() {
    }

    /**
     * Main entry point.
     */
    public static void main(final String[] args) throws Exception {
        final Args parsed = parseArgs(args);
        if (parsed.help) {
            printUsage();
            return;
        }

        if (parsed.outputFile == null) {
            throw new IllegalArgumentException("--output is required");
        }
        if (parsed.inputFile == null) {
            throw new IllegalArgumentException("--input is required");
        }
        if (!parsed.interleaved && parsed.inputFile2 == null) {
            throw new IllegalArgumentException("--input2 is required unless --interleaved is set");
        }
        if (parsed.interleaved && parsed.inputFile2 != null) {
            throw new IllegalArgumentException("--input2 is not allowed with --interleaved");
        }

        if (parsed.trainProp != null) {
            ClassifierFactory.setDataProp(parsed.trainProp, false);
        }

        final String gene = parsed.gene == null ? ClassifierFactory.RRNA_16S_GENE : parsed.gene;
        final ClassifierFactory factory = ClassifierFactory.getFactory(gene);
        final TrainingInfo trainingInfo = extractTrainingInfo(factory);
        final PairedNaiveBayesClassifier pairedClassifier = new PairedNaiveBayesClassifier(trainingInfo);

        final String[] ranks = useSpeciesRanks(gene)
                ? ClassificationResultFormatter.RANKS_WITHSPECIES
                : ClassificationResultFormatter.RANKS;

        try (PrintWriter outWriter = new PrintWriter(parsed.outputFile);
             PrintWriter shortSeqWriter = parsed.shortSeqFile == null ? null : new PrintWriter(parsed.shortSeqFile);
             PairReader pairReader = new PairReader(parsed.inputFile, parsed.inputFile2, parsed.interleaved)) {
            PairRecord pair;
            while ((pair = pairReader.nextPair()) != null) {
                try {
                    final ClassifierSequence leftSeq = new ClassifierSequence(pair.leftSeq);
                    final ClassifierSequence rightSeq = new ClassifierSequence(pair.rightSeq);
                    final ClassificationResult result = pairedClassifier.classifyPair(
                            pair.outputId,
                            leftSeq,
                            rightSeq,
                            parsed.minWords);
                    outWriter.print(formatResult(result, parsed.outputFormat, parsed.confidence, ranks));
                } catch (ShortSequenceException e) {
                    if (shortSeqWriter != null) {
                        shortSeqWriter.println(pair.outputId);
                    }
                }
            }
        }
    }

    /**
     * Render one classification result using stock formatter behavior.
     */
    private static String formatResult(
            final ClassificationResult result,
            final ClassificationResultFormatter.FORMAT outputFormat,
            final float confidence,
            final String[] ranks) {
        return ClassificationResultFormatter.getOutput(result, outputFormat, confidence, ranks);
    }

    /**
     * Use fungal ITS species-aware ranks when appropriate.
     */
    private static boolean useSpeciesRanks(final String gene) {
        return ClassifierFactory.FUNGALITS_warcup_GENE.equalsIgnoreCase(gene)
                || ClassifierFactory.FUNGALITS_unite_GENE.equalsIgnoreCase(gene);
    }

    /**
     * Reflection helper to reuse stock-loaded training objects without retraining.
     */
    private static TrainingInfo extractTrainingInfo(final ClassifierFactory factory) throws Exception {
        final Field trainingInfoField = ClassifierFactory.class.getDeclaredField("trainingInfo");
        trainingInfoField.setAccessible(true);
        final Object value = trainingInfoField.get(factory);
        return (TrainingInfo) value;
    }

    /**
     * Parse CLI arguments.
     */
    private static Args parseArgs(final String[] args) {
        final HashMap<String, String> valueOptions = new HashMap<String, String>();
        boolean interleaved = false;
        boolean help = false;

        int index = 0;
        while (index < args.length) {
            final String token = args[index];
            if ("--interleaved".equals(token)) {
                interleaved = true;
                index += 1;
                continue;
            }
            if ("--help".equals(token) || "-h".equals(token)) {
                help = true;
                index += 1;
                continue;
            }
            if (!token.startsWith("--")) {
                throw new IllegalArgumentException("Unexpected argument: " + token);
            }
            if (index + 1 >= args.length) {
                throw new IllegalArgumentException("Missing value for argument: " + token);
            }
            valueOptions.put(token, args[index + 1]);
            index += 2;
        }

        final Args parsed = new Args();
        parsed.help = help;
        parsed.interleaved = interleaved;
        parsed.inputFile = valueOptions.get("--input");
        parsed.inputFile2 = valueOptions.get("--input2");
        parsed.outputFile = valueOptions.get("--output");
        parsed.shortSeqFile = valueOptions.get("--shortseq-outfile");
        parsed.trainProp = valueOptions.get("--train-prop");
        parsed.gene = valueOptions.get("--gene");
        parsed.confidence = valueOptions.containsKey("--conf")
                ? Float.parseFloat(valueOptions.get("--conf"))
                : 0.8f;
        parsed.minWords = valueOptions.containsKey("--min-words")
                ? Integer.parseInt(valueOptions.get("--min-words"))
                : 5;

        final String formatText = valueOptions.containsKey("--format")
                ? valueOptions.get("--format")
                : "allrank";
        parsed.outputFormat = parseOutputFormat(formatText);

        if (ClassificationResultFormatter.FORMAT.biom.equals(parsed.outputFormat)) {
            throw new IllegalArgumentException("--format biom is not supported for paired TAV mode");
        }

        if (parsed.confidence < 0.0f || parsed.confidence > 1.0f) {
            throw new IllegalArgumentException("--conf must be in [0,1]");
        }
        if (parsed.minWords < 5) {
            throw new IllegalArgumentException("--min-words must be at least 5");
        }

        return parsed;
    }

    /**
     * Parse stock-compatible output format names.
     */
    private static ClassificationResultFormatter.FORMAT parseOutputFormat(final String formatText) {
        if ("allrank".equalsIgnoreCase(formatText)) {
            return ClassificationResultFormatter.FORMAT.allRank;
        }
        if ("fixrank".equalsIgnoreCase(formatText)) {
            return ClassificationResultFormatter.FORMAT.fixRank;
        }
        if ("filterbyconf".equalsIgnoreCase(formatText)) {
            return ClassificationResultFormatter.FORMAT.filterbyconf;
        }
        if ("db".equalsIgnoreCase(formatText)) {
            return ClassificationResultFormatter.FORMAT.dbformat;
        }
        if ("biom".equalsIgnoreCase(formatText)) {
            return ClassificationResultFormatter.FORMAT.biom;
        }
        throw new IllegalArgumentException("Unsupported --format value: " + formatText);
    }

    /**
     * Print command usage.
     */
    private static void printUsage() {
        System.out.println("USAGE: PairedClassifierMain [options]");
        System.out.println("  --input <file>               Required. R1 file or interleaved file");
        System.out.println("  --input2 <file>              Required unless --interleaved is set");
        System.out.println("  --interleaved                Interpret --input as interleaved paired file");
        System.out.println("  --output <file>              Required. TAV taxonomy output");
        System.out.println("  --gene <name>                16srrna|fungallsu|fungalits_warcup|fungalits_unite");
        System.out.println("  --train-prop <file>          Optional pretrained model properties override");
        System.out.println("  --format <name>              allrank|fixrank|filterbyconf|db (default allrank)");
        System.out.println("  --conf <float>               Confidence cutoff in [0,1] (default 0.8)");
        System.out.println("  --min-words <int>            Bootstrap min words, at least 5 (default 5)");
        System.out.println("  --shortseq-outfile <file>    Optional IDs of short/unclassifiable pairs");
        System.out.println("  --help                       Show this message");
    }

    /**
     * Parsed argument bag.
     */
    private static final class Args {
        private boolean help;
        private boolean interleaved;
        private String inputFile;
        private String inputFile2;
        private String outputFile;
        private String shortSeqFile;
        private String trainProp;
        private String gene;
        private float confidence;
        private int minWords;
        private ClassificationResultFormatter.FORMAT outputFormat;
    }

    /**
     * One synchronized paired sequence record.
     */
    private static final class PairRecord {
        private final Sequence leftSeq;
        private final Sequence rightSeq;
        private final String outputId;

        private PairRecord(final Sequence leftSeq, final Sequence rightSeq, final String outputId) {
            this.leftSeq = leftSeq;
            this.rightSeq = rightSeq;
            this.outputId = outputId;
        }
    }

    /**
     * Pair reader for two-file or interleaved input.
     */
    private static final class PairReader implements AutoCloseable {
        private final boolean interleaved;
        private final SeqReader readerLeft;
        private final SeqReader readerRight;

        private PairReader(
                final String inputFile,
                final String inputFile2,
                final boolean interleaved) throws IOException {
            this.interleaved = interleaved;
            this.readerLeft = new SequenceReader(new File(inputFile));
            if (interleaved) {
                this.readerRight = null;
            } else {
                this.readerRight = new SequenceReader(new File(inputFile2));
            }
        }

        /**
         * Read one paired record.
         */
        private PairRecord nextPair() throws IOException {
            final Sequence left;
            final Sequence right;

            if (interleaved) {
                left = readerLeft.readNextSequence();
                if (left == null) {
                    return null;
                }
                right = readerLeft.readNextSequence();
                if (right == null) {
                    throw new IllegalArgumentException("Interleaved input has odd number of records");
                }
            } else {
                left = readerLeft.readNextSequence();
                right = readerRight.readNextSequence();
                if (left == null && right == null) {
                    return null;
                }
                if (left == null || right == null) {
                    throw new IllegalArgumentException("Input files have different number of records");
                }
            }

            final String leftToken = firstToken(left.getSeqName());
            final String rightToken = firstToken(right.getSeqName());
            final String leftComparableId = comparableId(leftToken);
            final String rightComparableId = comparableId(rightToken);
            if (!leftComparableId.equals(rightComparableId)) {
                throw new IllegalArgumentException(
                        "Paired sequence ID mismatch: " + left.getSeqName() + " vs " + right.getSeqName());
            }
            return new PairRecord(left, right, leftToken);
        }

        /**
         * Close underlying readers.
         */
        @Override
        public void close() throws IOException {
            readerLeft.close();
            if (readerRight != null) {
                readerRight.close();
            }
        }

        /**
         * First whitespace-delimited token from a header.
         */
        private static String firstToken(final String rawHeader) {
            final String trimmed = rawHeader.trim();
            final int splitAt = trimmed.indexOf(' ');
            if (splitAt < 0) {
                return trimmed;
            }
            return trimmed.substring(0, splitAt);
        }

        /**
         * Comparable record ID for pair matching.
         */
        private static String comparableId(final String headerToken) {
            String normalized = headerToken;
            if (normalized.endsWith("/1") || normalized.endsWith("/2")) {
                normalized = normalized.substring(0, normalized.length() - 2);
            }
            final int semicolon = normalized.indexOf(';');
            if (semicolon >= 0) {
                normalized = normalized.substring(0, semicolon);
            }
            return normalized;
        }
    }
}
