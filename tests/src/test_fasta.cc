#include "fasta.h"
#include "gtest\gtest.h"

TEST(TestExtractRecords, ParseValidFile_FindAllRecords)
{
    auto extracted_records = mipfinder::fasta::extractRecords("valid.fasta");
    EXPECT_EQ(extracted_records.size(), 5);
}

TEST(TestExtractRecords, ParseFileWithGoodAndBadRecords_FindAllValidRecords)
{
    auto extracted_records = mipfinder::fasta::extractRecords("mixed.fasta");
    EXPECT_EQ(extracted_records.size(), 6);
}

TEST(TestExtractRecords, ParseFileWithNoRecords_FindNoValidRecords)
{
    auto extracted_records = mipfinder::fasta::extractRecords("empty.fasta");
    EXPECT_EQ(extracted_records.size(), 0);
}

TEST(TestExtractRecords, ParseRecordWithInvalidHeader_FindNoValidRecords)
{
    auto extracted_records = mipfinder::fasta::extractRecords("empty_header.fasta");
    EXPECT_EQ(extracted_records.size(), 0);
}

TEST(TestExtractRecords, ParseNonExistentFile_ThrowStdRuntimeError)
{
    EXPECT_THROW(mipfinder::fasta::extractRecords("does_not_exist.fasta"), std::runtime_error);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();
}