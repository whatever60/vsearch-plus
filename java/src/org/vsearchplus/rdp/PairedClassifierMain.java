package org.vsearchplus.rdp;

import edu.msu.cme.rdp.classifier.ClassificationResult;
import edu.msu.cme.rdp.classifier.HierarchyTree;
import edu.msu.cme.rdp.classifier.RankAssignment;
import edu.msu.cme.rdp.classifier.ShortSequenceException;
import edu.msu.cme.rdp.classifier.TrainingInfo;
import edu.msu.cme.rdp.classifier.io.ClassificationResultFormatter;
import edu.msu.cme.rdp.classifier.utils.ClassifierFactory;
import edu.msu.cme.rdp.classifier.utils.ClassifierSequence;
import edu.msu.cme.rdp.multicompare.MCSample;
import edu.msu.cme.rdp.multicompare.MCSamplePrintUtil;
import edu.msu.cme.rdp.multicompare.taxon.MCTaxon;
import edu.msu.cme.rdp.multicompare.visitors.DefaultPrintVisitor;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.taxatree.ConcretRoot;
import edu.msu.cme.rdp.taxatree.Taxon;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

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
        final PairedNaiveBayesClassifier pairedClassifier =
                new PairedNaiveBayesClassifier(trainingInfo, parsed.seed);
        final boolean hasCopyNumber = factory.getRoot().hasCopyNumberInfo();

        final String[] ranks = useSpeciesRanks(gene)
                ? ClassificationResultFormatter.RANKS_WITHSPECIES
                : ClassificationResultFormatter.RANKS;
        final boolean needsAggregateOutputs =
                (parsed.bootstrapOutFile != null) || (parsed.hierOutFile != null);
        final MCSample aggregateSample = needsAggregateOutputs ? new MCSample("paired_tav") : null;
        final List<MCSample> aggregateSamples = new ArrayList<MCSample>();
        if (aggregateSample != null) {
            aggregateSamples.add(aggregateSample);
        }
        final ConcretRoot<MCTaxon> hierRoot = (parsed.hierOutFile != null)
                ? new ConcretRoot<MCTaxon>(
                        new MCTaxon(factory.getRoot().getTaxid(), factory.getRoot().getName(), factory.getRoot().getRank()))
                : null;
        final HashMap<String, Long> seqCountMap = new HashMap<String, Long>();
        final HashMap<Taxon, Double> cachedCopyNumber = new HashMap<Taxon, Double>();

        try (PrintWriter outWriter = new PrintWriter(parsed.outputFile);
             PrintWriter shortSeqWriter = parsed.shortSeqFile == null ? null : new PrintWriter(parsed.shortSeqFile);
             PairReader pairReader = new PairReader(parsed.inputFile, parsed.inputFile2, parsed.interleaved)) {
            if (ClassificationResultFormatter.FORMAT.filterbyconf.equals(parsed.outputFormat)) {
                for (int i = 0; i < ranks.length; i++) {
                    outWriter.print("\t" + ranks[i]);
                }
                outWriter.println();
            }
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
                    if (aggregateSample != null) {
                        aggregateSample.addRankCount(result);
                    }
                    if (hierRoot != null) {
                        processClassificationResult(
                                result,
                                aggregateSample,
                                hierRoot,
                                parsed.confidence,
                                seqCountMap,
                                hasCopyNumber,
                                cachedCopyNumber);
                    }
                } catch (ShortSequenceException e) {
                    if (shortSeqWriter != null) {
                        shortSeqWriter.println(pair.outputId);
                    }
                }
            }
        }

        if (parsed.hierOutFile != null) {
            final File hierOutFile = new File(parsed.hierOutFile);
            try (PrintStream hierOut = new PrintStream(hierOutFile)) {
                final DefaultPrintVisitor printVisitor = new DefaultPrintVisitor(hierOut, aggregateSamples);
                hierRoot.topDownVisit(printVisitor);
            }

            if (hasCopyNumber) {
                final File parent = hierOutFile.getParentFile();
                final File copyNumberAdjusted = (parent == null)
                        ? new File("cnadjusted_" + hierOutFile.getName())
                        : new File(parent, "cnadjusted_" + hierOutFile.getName());
                try (PrintStream adjustedOut = new PrintStream(copyNumberAdjusted)) {
                    final DefaultPrintVisitor adjustedVisitor = new DefaultPrintVisitor(adjustedOut, aggregateSamples, true);
                    hierRoot.topDownVisit(adjustedVisitor);
                }
            }
        }

        if (parsed.bootstrapOutFile != null) {
            try (PrintStream bootstrapOut = new PrintStream(parsed.bootstrapOutFile)) {
                MCSamplePrintUtil.printBootstrapCountTable(bootstrapOut, aggregateSample);
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
        boolean queryFileFlag = false;

        int index = 0;
        while (index < args.length) {
            final String token = args[index];
            if ("--interleaved".equals(token)) {
                interleaved = true;
                index += 1;
                continue;
            }
            if ("--help".equals(token)) {
                help = true;
                index += 1;
                continue;
            }
            if ("--queryFile".equals(token) || "-q".equals(token)) {
                queryFileFlag = true;
                index += 1;
                continue;
            }

            final String normalized = normalizeOption(token);
            if (normalized == null) {
                throw new IllegalArgumentException("Unexpected argument: " + token);
            }
            if (index + 1 >= args.length) {
                throw new IllegalArgumentException("Missing value for argument: " + token);
            }
            valueOptions.put(normalized, args[index + 1]);
            index += 2;
        }

        final Args parsed = new Args();
        parsed.help = help;
        parsed.interleaved = interleaved;
        parsed.queryFile = queryFileFlag;
        parsed.inputFile = valueOptions.get("--input");
        parsed.inputFile2 = valueOptions.get("--input2");
        parsed.outputFile = valueOptions.get("--output");
        parsed.shortSeqFile = valueOptions.get("--shortseq-outfile");
        parsed.bootstrapOutFile = valueOptions.get("--bootstrap_outfile");
        parsed.hierOutFile = valueOptions.get("--hier_outfile");
        parsed.trainProp = valueOptions.get("--train-prop");
        parsed.gene = valueOptions.get("--gene");
        parsed.confidence = valueOptions.containsKey("--conf")
                ? Float.parseFloat(valueOptions.get("--conf"))
                : 0.8f;
        parsed.minWords = valueOptions.containsKey("--min-words")
                ? Integer.parseInt(valueOptions.get("--min-words"))
                : 5;
        parsed.seed = valueOptions.containsKey("--seed")
                ? Long.parseLong(valueOptions.get("--seed"))
                : 1L;

        final String formatText = valueOptions.containsKey("--format")
                ? valueOptions.get("--format")
                : "allrank";
        parsed.outputFormat = parseOutputFormat(formatText);

        if (ClassificationResultFormatter.FORMAT.biom.equals(parsed.outputFormat)) {
            throw new IllegalArgumentException("--format biom is not supported for paired TAV mode");
        }
        if (valueOptions.containsKey("--metadata")) {
            throw new IllegalArgumentException("--metadata is not supported for paired TAV mode");
        }
        if (valueOptions.containsKey("--biomFile")) {
            throw new IllegalArgumentException("--biomFile is not supported for paired TAV mode");
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
     * Normalize stock and extension option aliases to canonical keys.
     */
    private static String normalizeOption(final String token) {
        if ("--input".equals(token)) {
            return "--input";
        }
        if ("--input2".equals(token)) {
            return "--input2";
        }
        if ("--output".equals(token) || "--outputFile".equals(token) || "-o".equals(token)) {
            return "--output";
        }
        if ("--shortseq-outfile".equals(token) || "--shortseq_outfile".equals(token) || "-s".equals(token)) {
            return "--shortseq-outfile";
        }
        if ("--train-prop".equals(token) || "--train_propfile".equals(token) || "-t".equals(token)) {
            return "--train-prop";
        }
        if ("--gene".equals(token) || "-g".equals(token)) {
            return "--gene";
        }
        if ("--conf".equals(token) || "-c".equals(token)) {
            return "--conf";
        }
        if ("--format".equals(token) || "-f".equals(token)) {
            return "--format";
        }
        if ("--min-words".equals(token) || "--minWords".equals(token) || "-w".equals(token)) {
            return "--min-words";
        }
        if ("--seed".equals(token)) {
            return "--seed";
        }
        if ("--bootstrap_outfile".equals(token) || "-b".equals(token)) {
            return "--bootstrap_outfile";
        }
        if ("--hier_outfile".equals(token) || "-h".equals(token)) {
            return "--hier_outfile";
        }
        if ("--metadata".equals(token) || "-d".equals(token)) {
            return "--metadata";
        }
        if ("--biomFile".equals(token) || "-m".equals(token)) {
            return "--biomFile";
        }
        return null;
    }

    /**
     * Stock-style helper: find or create one taxon node under the aggregate output tree.
     */
    private static MCTaxon findOrCreateTaxon(
            final ConcretRoot<MCTaxon> root,
            final RankAssignment assignment,
            final int parentId,
            final boolean unclassified,
            final Map<String, Long> seqCountMap,
            final String lineage) {
        int taxid = assignment.getTaxid();
        if (unclassified) {
            taxid = Taxon.getUnclassifiedId(taxid);
        }

        MCTaxon taxon = root.getChildTaxon(taxid);
        if (taxon == null) {
            taxon = new MCTaxon(assignment.getTaxid(), assignment.getName(), assignment.getRank(), unclassified);
            root.addChild(taxon, parentId);

            Long val = seqCountMap.get(taxon.getRank());
            if (val == null) {
                val = 0L;
            }
            seqCountMap.put(taxon.getRank(), val + 1L);
            taxon.setLineage(lineage + taxon.getName() + ";" + taxon.getRank() + ";");
        }

        return taxon;
    }

    /**
     * Stock-style helper: resolve cached copy number for lowest confident assignment taxon.
     */
    private static double findCopyNumber(
            final RankAssignment assignment,
            final Taxon taxon,
            final Map<Taxon, Double> cache) {
        Double copyNumber = cache.get(taxon);
        if (copyNumber != null) {
            return copyNumber;
        }

        HierarchyTree current = assignment.getBestClass();
        while (current != null) {
            if (current.getName().equalsIgnoreCase(taxon.getName())
                    && current.getRank().equalsIgnoreCase(taxon.getRank())) {
                copyNumber = current.getCopyNumber();
                cache.put(taxon, copyNumber);
                return copyNumber;
            }
            current = current.getParent();
        }

        return 1.0;
    }

    /**
     * Stock-style helper: update aggregate hierarchy tree and counts for one classification result.
     */
    private static void processClassificationResult(
            final ClassificationResult result,
            final MCSample sample,
            final ConcretRoot<MCTaxon> root,
            final float confidence,
            final Map<String, Long> seqCountMap,
            final boolean hasCopyNumber,
            final Map<Taxon, Double> cachedCopyNumber) {
        RankAssignment lastAssignment = null;
        RankAssignment twoAgo = null;
        final StringBuilder lineage = new StringBuilder();
        MCTaxon taxon = null;
        MCTaxon countTaxon = null;
        final HashSet<MCTaxon> touchedTaxa = new HashSet<MCTaxon>();
        int parentId = root.getRootTaxid();
        final int count = sample.getDupCount(result.getSequence().getSeqName());

        for (RankAssignment assignment : result.getAssignments()) {
            boolean stop = false;
            if (assignment.getConfidence() < confidence) {
                parentId = root.getRootTaxid();
                if (twoAgo != null) {
                    parentId = twoAgo.getTaxid();
                }
                countTaxon = taxon;
                taxon = findOrCreateTaxon(root, lastAssignment, parentId, true, seqCountMap, lineage.toString());
                stop = true;
            } else {
                if (lastAssignment != null) {
                    parentId = lastAssignment.getTaxid();
                }
                taxon = findOrCreateTaxon(root, assignment, parentId, false, seqCountMap, lineage.toString());
                countTaxon = taxon;
            }

            touchedTaxa.add(taxon);
            twoAgo = lastAssignment;
            lastAssignment = assignment;

            if (stop) {
                break;
            }
            lineage.append(assignment.getName()).append(";").append(assignment.getRank()).append(";");
        }

        if (hasCopyNumber && (countTaxon != null) && (!result.getAssignments().isEmpty())) {
            final RankAssignment last = result.getAssignments().get(result.getAssignments().size() - 1);
            final double copyNumber = findCopyNumber(last, countTaxon, cachedCopyNumber);
            for (MCTaxon touched : touchedTaxa) {
                touched.incCount(sample, count, copyNumber);
            }
        } else {
            for (MCTaxon touched : touchedTaxa) {
                touched.incCount(sample, count);
            }
        }
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
        System.out.println("  -o,--output,--outputFile <file>");
        System.out.println("                                Required. TAV taxonomy output");
        System.out.println("  -g,--gene <name>             16srrna|fungallsu|fungalits_warcup|fungalits_unite");
        System.out.println("  -t,--train-prop,--train_propfile <file>");
        System.out.println("                                Optional pretrained model properties override");
        System.out.println("  -f,--format <name>           allrank|fixrank|filterbyconf|db (default allrank)");
        System.out.println("  -c,--conf <float>            Confidence cutoff in [0,1] (default 0.8)");
        System.out.println("  -w,--min-words,--minWords <int>");
        System.out.println("                                Bootstrap min words, at least 5 (default 5)");
        System.out.println("  --seed <int>                 Bootstrap random seed (default 1)");
        System.out.println("  -s,--shortseq-outfile,--shortseq_outfile <file>");
        System.out.println("                                Optional IDs of short/unclassifiable pairs");
        System.out.println("  -q,--queryFile               Legacy stock flag, accepted and ignored");
        System.out.println("  -b,--bootstrap_outfile <file>");
        System.out.println("                                Optional bootstrap summary output");
        System.out.println("  -h,--hier_outfile <file>");
        System.out.println("                                Optional hierarchy count output");
        System.out.println("  -d,--metadata                Parsed for parity, unsupported in paired TAV mode");
        System.out.println("  -m,--biomFile                Parsed for parity, unsupported in paired TAV mode");
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
        private String bootstrapOutFile;
        private String hierOutFile;
        private String trainProp;
        private String gene;
        private float confidence;
        private int minWords;
        private long seed;
        private ClassificationResultFormatter.FORMAT outputFormat;
        private boolean queryFile;
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
