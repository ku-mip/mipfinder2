#include "fasta.h"
#include "gtest\gtest.h"

TEST(TestExtractRecords, ParseValidFile_FindAllRecords)
{
    auto extracted_records = mipfinder::fasta::extractRecords("valid.fasta");
    EXPECT_EQ(extracted_records.size(), 5);
}

TEST(TestExtractRecords, ParseFileWithBadRecords_FindAllValidRecords)
{
    auto extracted_records = mipfinder::fasta::extractRecords("mixed.fasta");
    EXPECT_EQ(extracted_records.size(), 6);
}



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();
    std::cin.get();
}