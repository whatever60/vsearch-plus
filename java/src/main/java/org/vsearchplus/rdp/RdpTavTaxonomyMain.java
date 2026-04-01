package org.vsearchplus.rdp;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Java launcher for paired-end TAV taxonomy assignment.
 *
 * This bootstrap class keeps the Python launcher responsibilities in Java:
 * resolve the downloaded RDP runtime, compile the in-repo paired extension,
 * and launch the classifier main class with the correct classpath.
 */
public final class RdpTavTaxonomyMain {
  private static final String DEFAULT_RDP_ROOT =
      "data/rdp_classifier";
  private static final String DEFAULT_GENE = "16srrna";

  private RdpTavTaxonomyMain() {}

  /**
   * Main entry point.
   */
  public static void main(final String[] args) throws Exception {
    final LauncherArgs parsed = parseArgs(args);
    if (parsed.help) {
      printUsage();
      return;
    }

    final Path repoRoot = resolveRepoRoot();
    final Path classifierJar = resolveClassifierJar(parsed, repoRoot);
    final String trainProp = resolveEffectiveTrainProp(parsed, repoRoot);
    final List<String> runtimeJars = collectRuntimeJars(classifierJar);
    final Path buildClasses = compileJavaExtension(repoRoot, runtimeJars,
                                                   parsed.javacBin);

    final List<String> javaCommand =
        buildJavaCommand(parsed, buildClasses, runtimeJars, trainProp);
    final Process process = new ProcessBuilder(javaCommand)
                                .directory(repoRoot.toFile())
                                .inheritIO()
                                .start();
    final int exitCode = process.waitFor();
    if (exitCode != 0) {
      System.exit(exitCode);
    }

    System.out.println("RDP jar: " + classifierJar);
    if (trainProp != null) {
      System.out.println("Pretrained model property: " + trainProp);
    } else {
      System.out.println("Gene model: " + parsed.gene);
    }
    if (parsed.outputFile != null) {
      System.out.println("Output: " + parsed.outputFile);
    }
  }

  /**
   * Parse launcher-only arguments and retain child classifier arguments.
   */
  private static LauncherArgs parseArgs(final String[] args) {
    final LauncherArgs parsed = new LauncherArgs();

    int index = 0;
    while (index < args.length) {
      final String token = args[index];
      final String inlineValue = inlineValue(token);
      final String optionName = optionName(token);

      if ("--rdp-root".equals(optionName)) {
        parsed.rdpRoot = optionValue(args, index, token, inlineValue);
        index += inlineValue == null ? 2 : 1;
        continue;
      }
      if ("--rdp-jar".equals(optionName)) {
        parsed.rdpJar = optionValue(args, index, token, inlineValue);
        index += inlineValue == null ? 2 : 1;
        continue;
      }
      if ("--java-opt".equals(optionName)) {
        final String javaOpt = optionValue(args, index, token, inlineValue);
        parsed.javaOptions.addAll(splitShellWords(javaOpt));
        index += inlineValue == null ? 2 : 1;
        continue;
      }
      if ("--java-bin".equals(optionName)) {
        parsed.javaBin = optionValue(args, index, token, inlineValue);
        index += inlineValue == null ? 2 : 1;
        continue;
      }
      if ("--javac-bin".equals(optionName)) {
        parsed.javacBin = optionValue(args, index, token, inlineValue);
        index += inlineValue == null ? 2 : 1;
        continue;
      }

      parsed.childArgs.add(token);
      index += 1;
    }

    inspectChildArgs(parsed);
    return parsed;
  }

  /**
   * Inspect classifier arguments for help, gene, and explicit train props.
   */
  private static void inspectChildArgs(final LauncherArgs parsed) {
    int index = 0;
    while (index < parsed.childArgs.size()) {
      final String token = parsed.childArgs.get(index);
      final String optionName = optionName(token);
      final String inlineValue = inlineValue(token);

      if ("--help".equals(optionName)) {
        parsed.help = true;
        index += 1;
        continue;
      }
      if ("-t".equals(token) || "--train-prop".equals(optionName) ||
          "--train_propfile".equals(optionName)) {
        parsed.hasExplicitTrainProp = true;
        parsed.explicitTrainProp = childOptionValue(parsed.childArgs, index,
                                                    token, inlineValue);
        index += inlineValue == null ? 2 : 1;
        continue;
      }
      if ("-g".equals(token) || "--gene".equals(optionName)) {
        parsed.gene = childOptionValue(parsed.childArgs, index, token,
                                       inlineValue);
        index += inlineValue == null ? 2 : 1;
        continue;
      }
      if ("-o".equals(token) || "--output".equals(optionName) ||
          "--outputFile".equals(optionName)) {
        parsed.outputFile = childOptionValue(parsed.childArgs, index, token,
                                             inlineValue);
        index += inlineValue == null ? 2 : 1;
        continue;
      }

      index += 1;
    }
  }

  /**
   * Resolve the repository root from the launcher property or cwd.
   */
  private static Path resolveRepoRoot() {
    final String repoRoot = System.getProperty("vsearchplus.repo.root");
    if (repoRoot != null) {
      return Paths.get(repoRoot).toAbsolutePath().normalize();
    }
    return Paths.get("").toAbsolutePath().normalize();
  }

  /**
   * Resolve the stock RDP classifier jar path.
   */
  private static Path resolveClassifierJar(final LauncherArgs parsed,
                                           final Path repoRoot)
      throws IOException {
    if (parsed.rdpJar != null) {
      return resolvePath(repoRoot, parsed.rdpJar);
    }

    final Path manifestPath =
        resolvePath(repoRoot, parsed.rdpRoot).resolve("manifest.json");
    final Path classifierRoot =
        resolvePath(repoRoot, readManifestPath(manifestPath, "classifier"));

    try (Stream<Path> walk = Files.walk(classifierRoot)) {
      final List<Path> matches =
          walk.filter(path -> path.getFileName().toString().equals(
                           "classifier.jar"))
              .sorted()
              .collect(Collectors.toList());
      if (matches.isEmpty()) {
        throw new IllegalArgumentException("Could not find classifier.jar under "
                                           + classifierRoot);
      }
      return matches.get(0);
    }
  }

  /**
   * Resolve the effective train-prop path when one is available.
   */
  private static String resolveEffectiveTrainProp(final LauncherArgs parsed,
                                                  final Path repoRoot)
      throws IOException {
    if (parsed.hasExplicitTrainProp) {
      return parsed.explicitTrainProp;
    }

    final Path manifestPath =
        resolvePath(repoRoot, parsed.rdpRoot).resolve("manifest.json");
    final Path pretrainedRoot =
        resolvePath(repoRoot, readManifestPath(manifestPath, "pretrained_data"));
    final String relativeSuffix = String.join("/",
                                              "classifier",
                                              parsed.gene,
                                              "rRNAClassifier.properties");

    try (Stream<Path> walk = Files.walk(pretrainedRoot)) {
      final List<Path> matches =
          walk.filter(Files::isRegularFile)
              .filter(path -> path.toString().replace('\\', '/').endsWith(
                           relativeSuffix))
              .sorted()
              .collect(Collectors.toList());
      if (matches.isEmpty()) {
        return null;
      }
      return matches.get(0).toString();
    }
  }

  /**
   * Load one path value from the downloader manifest.
   */
  private static String readManifestPath(final Path manifestPath,
                                         final String key) throws IOException {
    boolean inPaths = false;
    for (String line : Files.readAllLines(manifestPath, StandardCharsets.UTF_8)) {
      final String trimmed = line.trim();
      if ("\"paths\": {".equals(trimmed)) {
        inPaths = true;
        continue;
      }
      if (inPaths && ("},".equals(trimmed) || "}".equals(trimmed))) {
        break;
      }
      final String prefix = "\"" + key + "\": ";
      if (inPaths && trimmed.startsWith(prefix)) {
        return unquoteJson(trimmed.substring(prefix.length()));
      }
    }
    throw new IllegalArgumentException("Could not find paths." + key + " in "
                                       + manifestPath);
  }

  /**
   * Resolve a possibly relative path from the repository root.
   */
  private static Path resolvePath(final Path repoRoot, final String rawPath) {
    final Path path = Paths.get(rawPath);
    if (path.isAbsolute()) {
      return path.normalize();
    }
    return repoRoot.resolve(path).normalize();
  }

  /**
   * Collect classifier and library jars needed for runtime compilation.
   */
  private static List<String> collectRuntimeJars(final Path classifierJar)
      throws IOException {
    final Path libDir = classifierJar.getParent().getParent().resolve("lib");
    final List<String> jars = new ArrayList<String>();
    jars.add(classifierJar.toString());

    if (Files.exists(libDir)) {
      try (Stream<Path> walk = Files.list(libDir)) {
        jars.addAll(
            walk.filter(path -> path.getFileName().toString().endsWith(".jar"))
                .sorted()
                .map(Path::toString)
                .collect(Collectors.toList()));
      }
    }
    return jars;
  }

  /**
   * Compile the in-repo paired classifier sources against the RDP jars.
   */
  private static Path compileJavaExtension(final Path repoRoot,
                                           final List<String> runtimeJars,
                                           final String javacBin)
      throws IOException, InterruptedException {
    final Path javaSrcDir =
        repoRoot.resolve("java").resolve("src").resolve("main").resolve(
            "java");
    final Path javaBuildDir =
        repoRoot.resolve("build").resolve("java").resolve("classes");
    Files.createDirectories(javaBuildDir);

    final List<String> sourceFiles;
    try (Stream<Path> walk = Files.walk(javaSrcDir)) {
      sourceFiles =
          walk.filter(path -> path.getFileName().toString().endsWith(".java"))
              .filter(path -> !path.getFileName().toString().equals(
                           "RdpTavTaxonomyMain.java"))
              .sorted()
              .map(Path::toString)
              .collect(Collectors.toList());
    }

    final List<String> compileCommand = new ArrayList<String>();
    compileCommand.add(javacBin);
    compileCommand.add("-cp");
    compileCommand.add(String.join(System.getProperty("path.separator"),
                                   runtimeJars));
    compileCommand.add("-d");
    compileCommand.add(javaBuildDir.toString());
    compileCommand.addAll(sourceFiles);

    final Process process = new ProcessBuilder(compileCommand)
                                .directory(repoRoot.toFile())
                                .inheritIO()
                                .start();
    final int exitCode = process.waitFor();
    if (exitCode != 0) {
      throw new IllegalStateException("javac failed with exit code " +
                                      exitCode);
    }
    return javaBuildDir;
  }

  /**
   * Build the child Java execution command.
   */
  private static List<String> buildJavaCommand(final LauncherArgs parsed,
                                               final Path buildClasses,
                                               final List<String> runtimeJars,
                                               final String trainProp) {
    final List<String> classpathEntries = new ArrayList<String>();
    classpathEntries.add(buildClasses.toString());
    classpathEntries.addAll(runtimeJars);

    final List<String> command = new ArrayList<String>();
    command.add(parsed.javaBin);
    command.addAll(parsed.javaOptions);
    command.add("-cp");
    command.add(String.join(System.getProperty("path.separator"),
                            classpathEntries));
    command.add("org.vsearchplus.rdp.PairedClassifierMain");
    command.addAll(parsed.childArgs);
    if ((!parsed.hasExplicitTrainProp) && (trainProp != null)) {
      command.add("--train-prop");
      command.add(trainProp);
    }
    return command;
  }

  /**
   * Shell-like splitting for repeated --java-opt values.
   */
  private static List<String> splitShellWords(final String text) {
    if (text.isEmpty()) {
      return Collections.emptyList();
    }

    final List<String> tokens = new ArrayList<String>();
    final StringBuilder current = new StringBuilder();
    boolean inSingle = false;
    boolean inDouble = false;
    boolean escaping = false;

    for (int index = 0; index < text.length(); index++) {
      final char ch = text.charAt(index);
      if (escaping) {
        current.append(ch);
        escaping = false;
        continue;
      }
      if ((ch == '\\') && (!inSingle)) {
        escaping = true;
        continue;
      }
      if ((ch == '\'') && (!inDouble)) {
        inSingle = !inSingle;
        continue;
      }
      if ((ch == '"') && (!inSingle)) {
        inDouble = !inDouble;
        continue;
      }
      if (Character.isWhitespace(ch) && (!inSingle) && (!inDouble)) {
        if (current.length() > 0) {
          tokens.add(current.toString());
          current.setLength(0);
        }
        continue;
      }
      current.append(ch);
    }

    if (escaping || inSingle || inDouble) {
      throw new IllegalArgumentException("Unterminated --java-opt value: " +
                                         text);
    }
    if (current.length() > 0) {
      tokens.add(current.toString());
    }
    return tokens;
  }

  /**
   * Extract the logical option name from a token.
   */
  private static String optionName(final String token) {
    final int equalsAt = token.indexOf('=');
    if (equalsAt < 0) {
      return token;
    }
    return token.substring(0, equalsAt);
  }

  /**
   * Extract the inline value from a token when present.
   */
  private static String inlineValue(final String token) {
    final int equalsAt = token.indexOf('=');
    if (equalsAt < 0) {
      return null;
    }
    return token.substring(equalsAt + 1);
  }

  /**
   * Resolve an option value for launcher-only arguments.
   */
  private static String optionValue(final String[] args, final int index,
                                    final String token,
                                    final String inlineValue) {
    if (inlineValue != null) {
      return inlineValue;
    }
    if (index + 1 >= args.length) {
      throw new IllegalArgumentException("Missing value for argument: " +
                                         token);
    }
    return args[index + 1];
  }

  /**
   * Resolve an option value from the forwarded child argument list.
   */
  private static String childOptionValue(final List<String> args,
                                         final int index, final String token,
                                         final String inlineValue) {
    if (inlineValue != null) {
      return inlineValue;
    }
    if (index + 1 >= args.size()) {
      throw new IllegalArgumentException("Missing value for argument: " +
                                         token);
    }
    return args.get(index + 1);
  }

  /**
   * Strip JSON quoting from one manifest line value.
   */
  private static String unquoteJson(final String rawValue) {
    String value = rawValue.trim();
    if (value.endsWith(",")) {
      value = value.substring(0, value.length() - 1);
    }
    if ((value.length() >= 2) && value.startsWith("\"") &&
        value.endsWith("\"")) {
      value = value.substring(1, value.length() - 1);
    }
    return value.replace("\\\\", "\\").replace("\\\"", "\"");
  }

  /**
   * Print combined launcher and classifier usage.
   */
  private static void printUsage() {
    System.out.println(
        "USAGE: vsearch-plus-rdp-tav [launcher options] [classification options]");
    System.out.println("Launcher options:");
    System.out.println(
        "  --rdp-root <dir>             RDP asset root containing manifest.json");
    System.out.println(
        "                                Default: " + DEFAULT_RDP_ROOT);
    System.out.println(
        "  --rdp-jar <file>             Override stock classifier.jar path");
    System.out.println(
        "  --java-opt <opt>             Extra JVM option for the classifier JVM");
    System.out.println(
        "                                Repeat to pass multiple options");
    System.out.println(
        "  --java-bin <path>            Java executable for classifier launch");
    System.out.println(
        "                                Default: current JVM executable");
    System.out.println(
        "  --javac-bin <path>           Javac executable for runtime compilation");
    System.out.println("                                Default: javac");
    System.out.println("  --help                       Show this message");
    System.out.println();
    System.out.println("Classification options:");
    System.out.println(
        "  --input <file>               Required. R1 file or interleaved file");
    System.out.println(
        "  --input2 <file>              Required unless --interleaved is set");
    System.out.println(
        "  --interleaved                Interpret --input as interleaved paired file");
    System.out.println("  -o,--output,--outputFile <file>");
    System.out.println(
        "                                Required. TAV taxonomy output");
    System.out.println(
        "  -g,--gene <name>             Default " + DEFAULT_GENE);
    System.out.println("  -t,--train-prop,--train_propfile <file>");
    System.out.println(
        "                                Optional pretrained model properties override");
    System.out.println(
        "  -f,--format <name>           allrank|fixrank|filterbyconf|db (default allrank)");
    System.out.println(
        "  -c,--conf <float>            Confidence cutoff in [0,1] (default 0.8)");
    System.out.println("  -w,--min-words,--minWords <int>");
    System.out.println(
        "                                Bootstrap min words, at least 5 (default 5)");
    System.out.println(
        "  --seed <int>                 Bootstrap random seed (default 1)");
    System.out.println("  -s,--shortseq-outfile,--shortseq_outfile <file>");
    System.out.println(
        "                                Optional IDs of short/unclassifiable pairs");
    System.out.println(
        "  -q,--queryFile               Legacy stock flag, accepted and ignored");
    System.out.println("  -b,--bootstrap_outfile <file>");
    System.out.println(
        "                                Optional bootstrap summary output");
    System.out.println("  -h,--hier_outfile <file>");
    System.out.println(
        "                                Optional hierarchy count output");
    System.out.println(
        "  -d,--metadata                Parsed for parity, unsupported in paired TAV mode");
    System.out.println(
        "  -m,--biomFile                Parsed for parity, unsupported in paired TAV mode");
  }

  /**
   * Parsed launcher argument bag.
   */
  private static final class LauncherArgs {
    private boolean help;
    private String rdpRoot = DEFAULT_RDP_ROOT;
    private String rdpJar;
    private String javaBin = defaultJavaBin();
    private String javacBin = "javac";
    private String gene = DEFAULT_GENE;
    private boolean hasExplicitTrainProp;
    private String explicitTrainProp;
    private String outputFile;
    private final List<String> javaOptions = new ArrayList<String>();
    private final List<String> childArgs = new ArrayList<String>();
  }

  /**
   * Resolve the current JVM binary path when available.
   */
  private static String defaultJavaBin() {
    final String javaHome = System.getProperty("java.home");
    if (javaHome == null) {
      return "java";
    }

    final boolean isWindows = System.getProperty("os.name")
                                  .toLowerCase(Locale.ROOT)
                                  .contains("win");
    final String binaryName = isWindows ? "java.exe" : "java";
    final Path candidate = Paths.get(javaHome).resolve("bin").resolve(
        binaryName);
    if (Files.isExecutable(candidate)) {
      return candidate.toString();
    }
    return "java";
  }
}
