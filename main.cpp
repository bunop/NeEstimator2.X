#include <cctype>
#include <cstdio>
#include <cstring>
#include <string>
#include <string_view>

extern "C" {
// C entry points implemented in Ne2x.c and reused by the C++ launcher.
int RunDirect(char misFilSuf[]);
int RunMultiFiles(char *mFileName, char mOpt);
int RunMultiCommon(char *mFileName);
int RunOption(char misFilSuf[], char LocSuf[], char BurSuf[], char hasOpt,
              char rem, char *FileOne, char *FileTwo);
}

// Anonymous namespace: everything inside this block is private to this file.
// It keeps helper functions from polluting the global symbol table and makes
// clear that they are internal implementation details of the C++ launcher.
namespace {

// Compare two strings without considering letter case.
// This is used so the C++ launcher preserves the legacy "rm" handling logic
// from the old dispatcher, where "rm", "RM", and "Rm" all behave alike.
bool equalsIgnoreCase(std::string_view left, std::string_view right) {
    if (left.size() != right.size()) {
        return false;
    }

    for (std::size_t index = 0; index < left.size(); ++index) {
        const unsigned char leftChar = static_cast<unsigned char>(left[index]);
        const unsigned char rightChar = static_cast<unsigned char>(right[index]);
        if (std::tolower(leftChar) != std::tolower(rightChar)) {
            return false;
        }
    }

    return true;
}

// Print the same error text as the original C dispatcher when the command line
// does not match the expected formats.
void printIllegalArgument() {
    std::printf("Illegal argument!\n");
}

// Reproduce the legacy summary line that reports how many input files were
// processed by the multi-file modes.
void printDataFileCount(int count) {
    if (count > 1) {
        std::printf("\n*** Number of data files = %d ***\n", count);
    }
}

// Remove a command prefix and keep only the file path or payload that follows
// it. Examples: i:path.txt -> path.txt, o:options.txt -> options.txt.
std::string tailAfterPrefix(std::string_view argument, std::size_t prefixLength) {
    if (argument.size() <= prefixLength) {
        return {};
    }

    return std::string(argument.substr(prefixLength));
}

}  // namespace

int main(int argc, char **argv) {
    // These suffixes are the same defaults used by the original C entry point.
    // They are passed unchanged into the existing C routines so the generated
    // output names and missing-data handling remain compatible.
    char misFilSuf[] = "NoDat.txt";
    char LocSuf[] = "Loc.txt";
    char BurSuf[] = "Bur.txt";

    // No argument means the program should enter the interactive/direct mode.
    if (argc < 2) {
        const int runs = RunDirect(misFilSuf);
        if (runs > 1) {
            std::printf("*** Number of runs = %d ***\n", runs);
        }
        return 0;
    }

    const std::string_view firstArgument = argv[1] != nullptr ? std::string_view(argv[1])
                                                              : std::string_view();
    const std::size_t firstLength = firstArgument.size();
    const char command = firstLength > 0 ? firstArgument.front() : '\0';

    // Only the legacy command families are accepted here:
    // i: for one info file, m: and m+: for multi-file runs, and c: for the
    // common-settings batch mode.
    if (firstLength <= 2 || (command != 'i' && command != 'm' && command != 'c')) {
        printIllegalArgument();
        return 0;
    }

    // i: and c: require an immediate colon after the command letter.
    if ((command == 'i' || command == 'c') && firstArgument[1] != ':') {
        printIllegalArgument();
        return 0;
    }

    char mOpt = 0;
    // The multi-file dispatcher has two forms:
    // m:  for the shorter legacy batch syntax
    // m+: for the extended syntax that enables extra parameters.
    if (command == 'm') {
        if (firstArgument[1] == '+') {
            if (firstLength == 3) {
                printIllegalArgument();
                return 0;
            }
            if (firstArgument[2] != ':') {
                printIllegalArgument();
                return 0;
            }
            mOpt = 1;
        } else {
            if (firstArgument[1] != ':') {
                printIllegalArgument();
                return 0;
            }
        }
    }

    std::string fileOne = tailAfterPrefix(firstArgument, mOpt == 1 ? 3U : 2U);

    // Multi-file mode: run the batch processor and optionally remove the
    // control file when the user asked for rm.
    if (command == 'm') {
        char rem = 0;
        if (argc > 2) {
            rem = equalsIgnoreCase(argv[2], "rm") ? 1 : 0;
        }

        const int dataFiles = RunMultiFiles(fileOne.data(), mOpt);
        printDataFileCount(dataFiles);

        if (rem == 1) {
            std::remove(fileOne.c_str());
        }
        return 0;
    }

    // Common-settings mode: use the shared configuration file and optionally
    // remove it afterward if rm was requested.
    if (command == 'c') {
        char rem = 0;
        if (argc > 2) {
            rem = equalsIgnoreCase(argv[2], "rm") ? 1 : 0;
        }

        const int dataFiles = RunMultiCommon(fileOne.data());
        printDataFileCount(dataFiles);

        if (rem == 1) {
            std::remove(fileOne.c_str());
        }
        return 0;
    }

    char hasOpt = 0;
    char rem = 0;
    std::string fileTwo;

    // Info-file mode: the first file is mandatory, the second file is optional
    // and starts with o:. The rm request can appear in either argument slot.
    if (argc > 2) {
        hasOpt = (argv[2][0] == 'o' && argv[2][1] == ':') ? 1 : 0;
        rem = (std::strcmp(argv[2], "rm") == 0) ? 1 : 0;
        if (hasOpt == 1 && argc > 3) {
            rem = (std::strcmp(argv[3], "rm") == 0) ? 1 : 0;
        }
    }

    if (hasOpt == 1) {
        fileTwo = tailAfterPrefix(std::string_view(argv[2]), 2U);
    }

    RunOption(misFilSuf, LocSuf, BurSuf, hasOpt, rem,
              fileOne.data(), fileTwo.data());
    return 0;
}
