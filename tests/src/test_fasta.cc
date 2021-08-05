#include "fasta.h"
#include "gtest\gtest.h"

TEST(TestExtractRecords, ParseValidFile_FindAllRecords)
{
    auto extracted_records = mipfinder::fasta::parse("valid.fasta");
    EXPECT_EQ(extracted_records.size(), 5);
}

TEST(TestExtractRecords, ParseFileWithGoodAndBadRecords_FindAllValidRecords)
{
    auto extracted_records = mipfinder::fasta::parse("mixed.fasta");
    EXPECT_EQ(extracted_records.size(), 6);
}

TEST(TestExtractRecords, ParseFileWithNoRecords_FindNoValidRecords)
{
    auto extracted_records = mipfinder::fasta::parse("empty.fasta");
    EXPECT_EQ(extracted_records.size(), 0);
}

TEST(TestExtractRecords, ParseRecordWithInvalidHeader_FindNoValidRecords)
{
    auto extracted_records = mipfinder::fasta::parse("empty_header.fasta");
    EXPECT_EQ(extracted_records.size(), 0);
}

TEST(TestExtractRecords, ParseNonExistentFile_ThrowStdRuntimeError)
{
    EXPECT_THROW(mipfinder::fasta::parse("does_not_exist.fasta"), std::runtime_error);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();
}